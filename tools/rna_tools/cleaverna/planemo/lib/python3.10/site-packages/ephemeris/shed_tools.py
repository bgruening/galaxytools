"""
A tool to automate installation of tool repositories from a Galaxy Tool Shed
into an instance of Galaxy.

Shed-tools has three commands: update, test and install.

Update simply updates all the tools in a Galaxy given connection details on the command line.

Test tests the specified tools in the Galaxy Instance.

Install allows installation of tools in multiple ways.
Galaxy instance details and the installed tools can be provided in one of three
ways:

1. In the YAML format via dedicated files (a sample can be found
   `here <https://github.com/galaxyproject/ansible-galaxy-tools/blob/master/files/tool_list.yaml.sample>`_).
2. On the command line as dedicated script options (see the usage help).
3. As a single composite parameter to the script. The parameter must be a
   single, YAML-formatted string with the keys corresponding to the keys
   available for use in the YAML formatted file (for example:
   `--yaml_tool "{'owner': 'kellrott', 'tool_shed_url':
   'https://testtoolshed.g2.bx.psu.edu', 'tool_panel_section_id':
   'peak_calling', 'name': 'synapse_interface'}"`).

Only one of the methods can be used with each invocation of the script but if
more than one are provided are provided, precedence will correspond to order
of the items in the list above.
When installing tools, Galaxy expects any `tool_panel_section_id` provided when
installing a tool to already exist in the configuration. If the section
does not exist, the tool will be installed outside any section. See
`shed_tool_conf.xml.sample` in this directory for a sample of such file. Before
running this script to install the tools, make sure to place such file into
Galaxy's configuration directory and set Galaxy configuration option
`tool_config_file` to include it.
"""

import datetime as dt
import json
import logging
import os
import re
import time
from collections import namedtuple
from collections.abc import Iterable
from concurrent.futures import (
    thread,
    ThreadPoolExecutor,
)
from typing import Optional

import requests
import yaml
from bioblend.galaxy.client import ConnectionError
from bioblend.galaxy.toolshed import ToolShedClient
from galaxy.tool_util.verify.interactor import (
    DictClientTestConfig,
    GalaxyInteractorApi,
    verify_tool,
)
from galaxy.util import unicodify
from typing_extensions import (
    NamedTuple,
    NotRequired,
    TypedDict,
)

from . import (
    get_galaxy_connection,
    load_yaml_file,
)
from .ephemeris_log import (
    disable_external_library_logging,
    setup_global_logger,
)
from .get_tool_list_from_galaxy import (
    GiToToolYaml,
    the_same_repository,
    tools_for_repository,
)
from .shed_tools_args import parser
from .shed_tools_methods import (
    complete_repo_information,
    flatten_repo_info,
    VALID_KEYS,
)

NON_TERMINAL_REPOSITORY_STATES = {
    "New",
    "Cloning",
    "Setting tool versions",
    "Installing repository dependencies",
    "Installing tool dependencies",
    "Loading proprietary datatypes",
}

log = logging.getLogger(__name__)


class InstallRepoDict(TypedDict):
    name: str
    owner: str
    changeset_revision: NotRequired[Optional[str]]
    tool_panel_section_id: NotRequired[Optional[str]]
    tool_panel_section_label: NotRequired[Optional[str]]
    tool_shed_url: NotRequired[str]
    revisions: NotRequired[list[str]]
    install_repository_dependencies: NotRequired[bool]
    install_resolver_dependencies: NotRequired[bool]
    install_tool_dependencies: NotRequired[bool]


class FilterResults(NamedTuple):
    already_installed_repos: list[InstallRepoDict]
    not_installed_repos: list[InstallRepoDict]


class InstallResults(NamedTuple):
    installed_repositories: list[InstallRepoDict]
    skipped_repositories: list[InstallRepoDict]
    errored_repositories: list[InstallRepoDict]


class InstallRepositoryManager:
    """Manages the installation of new repositories on a galaxy instance"""

    def __init__(self, galaxy_instance):
        """Initialize a new tool manager"""
        self.gi = galaxy_instance
        self.tool_shed_client = ToolShedClient(self.gi)

    def installed_repositories(self) -> list[InstallRepoDict]:
        """Get currently installed tools"""
        return GiToToolYaml(
            gi=self.gi,
            skip_tool_panel_section_name=False,
            get_data_managers=True,
            get_all_tools=True,
        ).tool_list.get("tools")

    def filter_installed_repos(self, repos: Iterable[InstallRepoDict], check_revision: bool = True) -> FilterResults:
        """This filters a list of repositories"""
        not_installed_repos: list[InstallRepoDict] = []
        already_installed_repos: list[InstallRepoDict] = []
        if check_revision:
            # If we want to check if revisions are equal, flatten the list,
            # so each repository - revision combination has its own entry
            installed_repos = flatten_repo_info(self.installed_repositories())
        else:
            # If we do not care about revision equality, do not do the flatten
            # action to limit the number of comparisons.
            installed_repos = self.installed_repositories()

        for repo in repos:
            for installed_repo in installed_repos:
                if the_same_repository(installed_repo, repo, check_revision):
                    already_installed_repos.append(repo)
                    break
            else:  # This executes when the for loop completes and no match has been found.
                not_installed_repos.append(repo)
        return FilterResults(
            already_installed_repos=already_installed_repos,
            not_installed_repos=not_installed_repos,
        )

    def install_repositories(
        self,
        repositories: list[InstallRepoDict],
        log=log,
        force_latest_revision: bool = False,
        default_toolshed: str = "https://toolshed.g2.bx.psu.edu/",
        default_install_tool_dependencies: bool = False,
        default_install_resolver_dependencies: bool = True,
        default_install_repository_dependencies: bool = True,
    ):
        """Install a list of tools on the current galaxy"""
        installation_start = dt.datetime.now()
        installed_repositories: list[InstallRepoDict] = []
        skipped_repositories: list[InstallRepoDict] = []
        errored_repositories: list[InstallRepoDict] = []
        counter = 0

        # Check repos for invalid keys
        for repo in repositories:
            for key in repo.keys():
                if key not in VALID_KEYS and key != "revisions":
                    if log:
                        log.warning(f"'{key}' not a valid key. Will be skipped during parsing")

        # Start by flattening the repo list per revision
        flattened_repos = flatten_repo_info(repositories)
        total_num_repositories = len(flattened_repos)

        # Complete the repo information, and make sure each repository has a revision
        repository_list: list[InstallRepoDict] = []
        for repository in flattened_repos:
            start = dt.datetime.now()
            try:
                complete_repo = complete_repo_information(
                    repository,
                    default_toolshed_url=default_toolshed,
                    default_install_tool_dependencies=default_install_tool_dependencies,
                    default_install_resolver_dependencies=default_install_resolver_dependencies,
                    default_install_repository_dependencies=default_install_repository_dependencies,
                    force_latest_revision=force_latest_revision,
                )
                repository_list.append(complete_repo)
            except Exception as e:
                # We'll run through the loop come whatever may, we log the errored repositories at the end anyway.
                if log:
                    log_repository_install_error(repository, start, unicodify(e), log)
                errored_repositories.append(repository)

        # Filter out already installed repos
        filtered_repos = self.filter_installed_repos(repository_list)

        for skipped_repo in filtered_repos.already_installed_repos:
            counter += 1
            if log:
                log_repository_install_skip(skipped_repo, counter, total_num_repositories, log)
            skipped_repositories.append(skipped_repo)

        # Install repos
        for repository in filtered_repos.not_installed_repos:
            counter += 1
            if log:
                log_repository_install_start(
                    repository,
                    counter=counter,
                    installation_start=installation_start,
                    log=log,
                    total_num_repositories=total_num_repositories,
                )
            result = self.install_repository_revision(repository, log)
            if result == "error":
                errored_repositories.append(repository)
            elif result == "skipped":
                skipped_repositories.append(repository)
            elif result == "installed":
                installed_repositories.append(repository)

        # Log results
        if log:
            log.info(
                "Installed repositories ({}): {}".format(
                    len(installed_repositories),
                    [(t["name"], t.get("changeset_revision")) for t in installed_repositories],
                )
            )
            log.info(
                "Skipped repositories ({}): {}".format(
                    len(skipped_repositories),
                    [(t["name"], t.get("changeset_revision")) for t in skipped_repositories],
                )
            )
            log.info(
                "Errored repositories ({}): {}".format(
                    len(errored_repositories),
                    [(t["name"], t.get("changeset_revision", "")) for t in errored_repositories],
                )
            )
            log.info("All repositories have been installed.")
            log.info(f"Total run time: {dt.datetime.now() - installation_start}")
        return InstallResults(
            installed_repositories=installed_repositories,
            skipped_repositories=skipped_repositories,
            errored_repositories=errored_repositories,
        )

    def update_repositories(self, repositories=None, log=log, **kwargs):
        if not repositories:  # Repositories None or empty list
            repositories = self.installed_repositories()
        else:
            filtered_repos = self.filter_installed_repos(repositories, check_revision=False)
            if filtered_repos.not_installed_repos:
                if log:
                    log.warning(
                        f"The following tools are not installed and will not be upgraded: {filtered_repos.not_installed_repos}"
                    )
            repositories = filtered_repos.already_installed_repos
        return self.install_repositories(repositories, force_latest_revision=True, log=log, **kwargs)

    def test_tools(
        self,
        test_json,
        repositories=None,
        log=log,
        test_user_api_key=None,
        test_user="ephemeris@galaxyproject.org",
        test_history_name=None,
        parallel_tests=1,
        test_all_versions=False,
        client_test_config_path=None,
    ):
        """Run tool tests for all tools in each repository in supplied tool list or ``self.installed_repositories()``."""
        tool_test_start = dt.datetime.now()
        tests_passed = []
        test_exceptions = []

        if not repositories:  # If repositories is None or empty list
            # Consider a variant of this that doesn't even consume a tool list YAML? target
            # something like installed_repository_revisions(self.gi)
            repositories = self.installed_repositories()

        target_repositories = flatten_repo_info(repositories)

        installed_tools = []
        for target_repository in target_repositories:
            repo_tools = tools_for_repository(self.gi, target_repository, all_tools=test_all_versions)
            installed_tools.extend(repo_tools)

        all_test_results = []
        galaxy_interactor = self._get_interactor(test_user, test_user_api_key)
        if client_test_config_path is not None:
            with open(client_test_config_path) as f:
                client_test_config_dict = yaml.full_load(f)
            client_test_config = DictClientTestConfig(client_test_config_dict.get("tools"))
        else:
            client_test_config = None

        if test_history_name:
            for history in self.gi.histories.get_histories(name=test_history_name, deleted=False):
                test_history = history["id"]
                log.debug(
                    "Using existing history with id '%s', last updated: %s",
                    test_history,
                    history["update_time"],
                )
                break
            else:
                test_history = galaxy_interactor.new_history(history_name=test_history_name)
        else:
            test_history = galaxy_interactor.new_history()

        with ThreadPoolExecutor(max_workers=parallel_tests) as executor:
            try:
                for tool in installed_tools:
                    self._test_tool(
                        executor=executor,
                        tool=tool,
                        galaxy_interactor=galaxy_interactor,
                        test_history=test_history,
                        log=log,
                        tool_test_results=all_test_results,
                        tests_passed=tests_passed,
                        test_exceptions=test_exceptions,
                        client_test_config=client_test_config,
                    )
            finally:
                # Always write report, even if test was cancelled.
                try:
                    executor.shutdown(wait=True)
                except KeyboardInterrupt:
                    executor._threads.clear()
                    thread._threads_queues.clear()
                n_passed = len(tests_passed)
                n_failed = len(test_exceptions)
                report_obj = {
                    "version": "0.1",
                    "suitename": f"Ephemeris tool tests targeting {self.gi.base_url}",
                    "results": {
                        "total": n_passed + n_failed,
                        "errors": n_failed,
                        "failures": 0,
                        "skips": 0,
                    },
                    "tests": sorted(all_test_results, key=lambda el: el["id"]),
                }
                with open(test_json, "w") as f:
                    json.dump(report_obj, f)
                if log:
                    log.info("Report written to '%s'", os.path.abspath(test_json))
                    log.info(f"Passed tool tests ({n_passed}): {[t for t in tests_passed]}")
                    log.info(f"Failed tool tests ({n_failed}): {[t[0] for t in test_exceptions]}")
                    log.info(f"Total tool test time: {dt.datetime.now() - tool_test_start}")

    def _get_interactor(self, test_user, test_user_api_key):
        if test_user_api_key is None:
            whoami = self.gi.make_get_request(self.gi.url + "/whoami").json()
            if whoami is not None:
                test_user_api_key = self.gi.key
        galaxy_interactor_kwds = {
            "galaxy_url": re.sub("/api", "", self.gi.url),
            "master_api_key": self.gi.key,
            "api_key": test_user_api_key,  # TODO
            "keep_outputs_dir": "",
        }
        if test_user_api_key is None:
            galaxy_interactor_kwds["test_user"] = test_user
        galaxy_interactor = GalaxyInteractorApi(**galaxy_interactor_kwds)
        return galaxy_interactor

    @staticmethod
    def _test_tool(
        executor,
        tool,
        galaxy_interactor,
        tool_test_results,
        tests_passed,
        test_exceptions,
        log,
        test_history=None,
        client_test_config=None,
    ):
        if test_history is None:
            test_history = galaxy_interactor.new_history()
        tool_id = tool["id"]
        tool_version = tool["version"]
        # If given a tool_id with a version suffix, strip it off so we can treat tool_version
        # correctly at least in client_test_config.
        if tool_version and tool_id.endswith("/" + tool_version):
            tool_id = tool_id[: -len("/" + tool_version)]

        label_base = tool_id
        if tool_version:
            label_base += "/" + str(tool_version)
        try:
            tool_test_dicts = galaxy_interactor.get_tool_tests(tool_id, tool_version=tool_version)
        except Exception as e:
            if log:
                log.warning(
                    "Fetching test definition for tool '%s' failed",
                    label_base,
                    exc_info=True,
                )
            test_exceptions.append((label_base, e))
            Results = namedtuple("Results", ["tool_test_results", "tests_passed", "test_exceptions"])
            return Results(
                tool_test_results=tool_test_results,
                tests_passed=tests_passed,
                test_exceptions=test_exceptions,
            )
        test_indices = list(range(len(tool_test_dicts)))

        for test_index in test_indices:
            test_id = label_base + "-" + str(test_index)

            def run_test(index, test_id):
                def register(job_data):
                    tool_test_results.append(
                        {
                            "id": test_id,
                            "has_data": True,
                            "data": job_data,
                        }
                    )

                try:
                    if log:
                        log.info("Executing test '%s'", test_id)
                    verify_tool(
                        tool_id,
                        galaxy_interactor,
                        test_index=index,
                        tool_version=tool_version,
                        register_job_data=register,
                        quiet=True,
                        test_history=test_history,
                        client_test_config=client_test_config,
                    )
                    tests_passed.append(test_id)
                    if log:
                        log.info("Test '%s' passed", test_id)
                except Exception as e:
                    if log:
                        log.warning("Test '%s' failed", test_id, exc_info=True)
                    test_exceptions.append((test_id, e))

            executor.submit(run_test, test_index, test_id)

    def install_repository_revision(self, repository: InstallRepoDict, log):
        default_err_msg = "All repositories that you are attempting to install " "have been previously installed."
        start = dt.datetime.now()
        try:
            response = self.tool_shed_client.install_repository_revision(
                tool_shed_url=repository["tool_shed_url"],
                name=repository["name"],
                owner=repository["owner"],
                changeset_revision=repository["changeset_revision"],
                install_tool_dependencies=repository["install_tool_dependencies"],
                install_repository_dependencies=repository["install_repository_dependencies"],
                install_resolver_dependencies=repository["install_resolver_dependencies"],
                new_tool_panel_section_label=repository.get("tool_panel_section_label"),
                tool_panel_section_id=repository.get("tool_panel_section_id"),
            )
            if isinstance(response, dict) and response.get("status", None) == "ok":
                # This rare case happens if a repository is already installed but
                # was not recognised as such in the above check. In such a
                # case the return value looks like this:
                # {u'status': u'ok', u'message': u'No repositories were
                #  installed, possibly because the selected repository has
                #  already been installed.'}
                if log:
                    log.debug("\tRepository {} is already installed.".format(repository["name"]))
            if log:
                log_repository_install_success(repository=repository, start=start, log=log)
            return "installed"
        except (ConnectionError, requests.exceptions.ConnectionError) as e:
            if default_err_msg in unicodify(e):
                # THIS SHOULD NOT HAPPEN DUE TO THE CHECKS EARLIER
                if log:
                    log.debug(
                        "\tRepository {} already installed (at revision {})".format(
                            repository["name"], repository["changeset_revision"]
                        )
                    )
                return "skipped"
            elif "504" in unicodify(e) or "Connection aborted" in unicodify(e):
                if log:
                    log.debug(
                        "Timeout during install of %s, extending wait to 1h",
                        repository["name"],
                    )
                success = self.wait_for_install(repository=repository, log=log, timeout=3600)
                if success:
                    if log:
                        log_repository_install_success(repository=repository, start=start, log=log)
                    return "installed"
                else:
                    if log:
                        log_repository_install_error(
                            repository=repository,
                            start=start,
                            msg=getattr(e, "body", unicodify(e)),
                            log=log,
                        )
                    return "error"
            else:
                if log:
                    log_repository_install_error(
                        repository=repository,
                        start=start,
                        msg=getattr(e, "body", unicodify(e)),
                        log=log,
                    )
                return "error"

    def wait_for_install(self, repository, log=log, timeout=3600):
        """
        If nginx times out, we look into the list of installed repositories
        and try to determine if a repository of the same namer/owner is still installing.
        Returns True if install finished successfully,
        returns False when timeout is exceeded or installation has failed.
        """
        # We request a repository revision, but Galaxy may decide to install the next downloable revision.
        # This ensures we have a revision to track, and if not, finds the revision that is actually being installed
        name = repository["name"]
        owner = repository["owner"]
        changeset_revision = repository["changeset_revision"]
        installed_repos = self.tool_shed_client.get_repositories()
        filtered_repos = [r for r in installed_repos if r["name"] == name and r["owner"] == owner]
        assert filtered_repos, f"Repository '{name}' from owner '{owner}' not in list of repositories."
        # Check if exact repository revision in filtered_repos
        installing_repo_id = None
        for repo in filtered_repos:
            if repo["changeset_revision"] == changeset_revision:
                installing_repo_id = repo["id"]
                break
        else:
            # Galaxy may have decided to install a newer repository revision. We now try to guess which repository that is.
            non_terminal = [r for r in filtered_repos if r["status"] in NON_TERMINAL_REPOSITORY_STATES]
            if len(non_terminal) == 1:
                # Unambiguous, we wait for this repo
                installing_repo_id = non_terminal[0]["id"]
            elif len(filtered_repos) == 1:
                installing_repo_id = filtered_repos[0]["id"]
            else:
                # We may have a repo that is permanently in a non-terminal state (e.g because of restart during installation).
                # Raise an exception and continue with the remaining repos.
                msg = "Could not track repository for name '%s', owner '%s', revision '%s'. "
                msg += "Please uninstall all non-terminal repositories and ensure revision '%s' is installable."
                raise AssertionError(msg % (name, owner, changeset_revision, changeset_revision))
        start = dt.datetime.now()
        while (dt.datetime.now() - start) < dt.timedelta(seconds=timeout):
            try:
                installed_repo = self.tool_shed_client.show_repository(installing_repo_id)
                status = installed_repo["status"]
                if status == "Installed":
                    return True
                elif status == "Error":
                    return False
                elif status in NON_TERMINAL_REPOSITORY_STATES:
                    time.sleep(10)
                else:
                    raise AssertionError(f"Repository name '{name}', owner '{owner}' in unknown status '{status}'")
            except ConnectionError as e:
                if log:
                    log.warning("Failed to get repositories list: %s", unicodify(e))
                time.sleep(10)
        return False


def log_repository_install_error(repository, start, msg, log):
    """
    Log failed repository installations. Return a dictionary with information
    """
    end = dt.datetime.now()
    log.error(
        "\t* Error installing a repository (after %s seconds)! Name: %s," "owner: %s, " "revision: %s, error: %s",
        str(end - start),
        repository.get("name", ""),
        repository.get("owner", ""),
        repository.get("changeset_revision", ""),
        msg,
    )


def log_repository_install_success(repository, start, log):
    """
    Log successful repository installation.
    Repositories that finish in error still count as successful installs currently.
    """
    end = dt.datetime.now()
    log.debug(
        "\trepository {} installed successfully (in {}) at revision {}".format(
            repository["name"], str(end - start), repository["changeset_revision"]
        )
    )


def log_repository_install_skip(repository, counter, total_num_repositories, log):
    log.debug(
        "({}/{}) repository {} already installed at revision {}. Skipping.".format(
            counter,
            total_num_repositories,
            repository["name"],
            repository["changeset_revision"],
        )
    )


def log_repository_install_start(
    repository: InstallRepoDict,
    counter,
    total_num_repositories,
    installation_start,
    log,
):
    log.debug(
        '({}/{}) Installing repository {} from {} to section "{}" at revision {} (TRT: {})'.format(
            counter,
            total_num_repositories,
            repository["name"],
            repository["owner"],
            repository.get("tool_panel_section_id") or repository.get("tool_panel_section_label"),
            repository.get("changeset_revision"),
            dt.datetime.now() - installation_start,
        )
    )


def args_to_repos(args) -> list[InstallRepoDict]:
    if args.tool_list_file:
        tool_list = load_yaml_file(args.tool_list_file)
        repos = tool_list["tools"]
    elif args.tool_yaml:
        repos = [yaml.safe_load(args.tool_yaml)]
    elif args.name and args.owner:
        repo = dict(
            owner=args.owner,
            name=args.name,
            tool_panel_section_id=args.tool_panel_section_id,
            tool_panel_section_label=args.tool_panel_section_label,
            revisions=args.revisions,
        )
        if args.tool_shed_url:
            repo["tool_shed_url"] = args.tool_shed_url
        repos = [repo]
    else:
        repos = []
    return repos


def main(argv=None):
    disable_external_library_logging()
    args = parser().parse_args(argv)
    log = setup_global_logger(name=__name__, log_file=args.log_file)
    gi = get_galaxy_connection(args, file=args.tool_list_file, log=log, login_required=True)
    install_repository_manager = InstallRepositoryManager(gi)

    repos = args_to_repos(args)

    if args.tool_list_file:
        tool_list = load_yaml_file(args.tool_list_file)
    else:
        tool_list = dict()

    # Get some of the other installation arguments
    kwargs = dict(
        default_install_tool_dependencies=tool_list.get("install_tool_dependencies")
        or getattr(args, "install_tool_dependencies", False),
        default_install_repository_dependencies=tool_list.get("install_repository_dependencies")
        or getattr(args, "install_repository_dependencies", False),
        default_install_resolver_dependencies=tool_list.get("install_resolver_dependencies")
        or getattr(args, "install_resolver_dependencies", False),
    )

    # Start installing/updating and store the results in install_results.
    # Or do testing if the action is `test`
    install_results = None
    if args.action == "update":
        install_results = install_repository_manager.update_repositories(repositories=repos, log=log, **kwargs)
    elif args.action == "install":
        install_results = install_repository_manager.install_repositories(
            repos, log=log, force_latest_revision=args.force_latest_revision, **kwargs
        )
    elif args.action == "test":
        install_repository_manager.test_tools(
            test_json=args.test_json,
            repositories=repos,
            log=log,
            test_user_api_key=args.test_user_api_key,
            test_user=args.test_user,
            test_history_name=args.test_history_name,
            parallel_tests=args.parallel_tests,
            test_all_versions=args.test_all_versions,
            client_test_config_path=args.client_test_config,
        )
    else:
        raise NotImplementedError("This point in the code should not be reached. Please contact the developers.")

    # Run tests on the install results if required.
    if install_results and args.test or args.test_existing:
        to_be_tested_repositories = install_results.installed_repositories
        if args.test_existing:
            to_be_tested_repositories.extend(install_results.skipped_repositories)
        if to_be_tested_repositories:
            install_repository_manager.test_tools(
                test_json=args.test_json,
                repositories=to_be_tested_repositories,
                log=log,
                test_user_api_key=args.test_user_api_key,
                test_user=args.test_user,
                parallel_tests=args.parallel_tests,
                client_test_config_path=args.client_test_config,
            )


if __name__ == "__main__":
    main()
