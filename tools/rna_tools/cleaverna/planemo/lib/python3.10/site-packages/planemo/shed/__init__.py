"""Abstractions for shed related interactions used by the rest of planemo."""

import contextlib
import copy
import fnmatch
import hashlib
import json
import os
import re
import shutil
import sys
import tarfile
from tempfile import mkstemp
from typing import (
    NamedTuple,
    TYPE_CHECKING,
)

import bioblend
import yaml
from bioblend.toolshed import ToolShedInstance
from galaxy.util import (
    odict,
    unicodify,
)

from planemo import (
    git,
    glob,
    templates,
)
from planemo.io import (
    can_write_to_path,
    coalesce_return_codes,
    error,
    find_matching_directories,
    info,
    shell,
    temp_directory,
    warn,
)
from planemo.shed2tap.base import BasePackage
from planemo.tools import yield_tool_sources
from .diff import diff_and_remove
from .interface import (
    api_exception_to_message,
    download_tar,
    find_category_ids,
    find_repository,
    latest_installable_revision,
    tool_shed_instance,
    username,
)

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext


SHED_CONFIG_NAME = ".shed.yml"
DOCKSTORE_REGISTRY_CONF = ".dockstore.yml"
REPO_METADATA_FILES = (SHED_CONFIG_NAME, DOCKSTORE_REGISTRY_CONF)
REPO_DEPENDENCIES_CONFIG_NAME = "repository_dependencies.xml"
TOOL_DEPENDENCIES_CONFIG_NAME = "tool_dependencies.xml"

NO_REPOSITORIES_MESSAGE = "Could not find any .shed.yml files or a --name to describe the target repository."
NAME_INVALID_MESSAGE = "Cannot use --name argument when multiple directories in target contain .shed.yml files."
NAME_REQUIRED_MESSAGE = "No repository name discovered but one is required."
CONFLICTING_NAMES_MESSAGE = "The supplied name argument --name conflicts with value discovered in .shed.yml."
PARSING_PROBLEM = "Problem parsing file .shed.yml in directory %s, skipping repository. Message: [%s]."
AUTO_REPO_CONFLICT_MESSAGE = "Cannot specify both auto_tool_repositories and repositories in .shed.yml at this time."
AUTO_NAME_CONFLICT_MESSAGE = (
    "Cannot specify both auto_tool_repositories and in .shed.yml and --name on the command-line."
)
REALIZAION_PROBLEMS_MESSAGE = "Problem encountered executing action for one or more repositories."
INCORRECT_OWNER_MESSAGE = (
    "Attempting to create a repository with configured owner [%s] that does not match API user [%s]."
)
PROBLEM_PROCESSING_REPOSITORY_MESSAGE = "Problem processing repositories, exiting."

# Planemo generated or consumed files that do not need to be uploaded to the
# tool shed.
PLANEMO_FILES = [
    "shed_upload*.tar.gz",
    "shed_download*.tar.gz",
    "tool_test_output.*",
    ".travis",
    ".travis.yml",
    ".shed.yml",
    "*~",
    "#*#",
]
SHED_SHORT_NAMES = {
    "toolshed": "https://toolshed.g2.bx.psu.edu/",
    "testtoolshed": "https://testtoolshed.g2.bx.psu.edu/",
    "local": "http://localhost:9009/",
}
SHED_LABELS = {
    "toolshed": "main Tool Shed",
    "testtoolshed": "test Tool Shed",
    "local": "local Tool Shed",
}
REPO_TYPE_UNRESTRICTED = "unrestricted"
REPO_TYPE_TOOL_DEP = "tool_dependency_definition"
REPO_TYPE_SUITE = "repository_suite_definition"

# TODO: sync this with tool shed impl someday
VALID_REPOSITORYNAME_RE = re.compile(r"^[a-z0-9\_]+$")
VALID_PUBLICNAME_RE = re.compile(r"^[a-z0-9._\-]+$")


# Generate with python scripts/categories.py
CURRENT_CATEGORIES = [
    "Assembly",
    "Astronomy",
    "ChIP-seq",
    "Climate Analysis",
    "CLIP-seq",
    "Combinatorial Selections",
    "Computational chemistry",
    "Constructive Solid Geometry",
    "Convert Formats",
    "Data Export",
    "Data Managers",
    "Data Source",
    "Ecology",
    "Entomology",
    "Epigenetics",
    "Fasta Manipulation",
    "Fastq Manipulation",
    "Flow Cytometry Analysis",
    "Genome annotation",
    "Genome editing",
    "Genome-Wide Association Study",
    "Genomic Interval Operations",
    "Geo Science",
    "GIS",
    "Graphics",
    "Imaging",
    "Interactive Tools",
    "Machine Learning",
    "Materials science",
    "Metabolomics",
    "Metagenomics",
    "Micro-array Analysis",
    "Molecular Dynamics",
    "Muon spectroscopy",
    "Nanopore",
    "Natural Language Processing",
    "Next Gen Mappers",
    "Ontology Manipulation",
    "Phylogenetics",
    "Proteomics",
    "RNA",
    "SAM",
    "Sequence Analysis",
    "Single Cell",
    "Spatial Omics",
    "Statistics",
    "Structural Materials Analysis",
    "Synthetic Biology",
    "Systems Biology",
    "Text Manipulation",
    "Tool Dependency Packages",
    "Tool Generators",
    "Transcriptomics",
    "Variant Analysis",
    "Visualization",
    "Web Services",
]
# http://stackoverflow.com/questions/7676255/find-and-replace-urls-in-a-block-of-te
HTTP_REGEX_PATTERN = re.compile(
    r"""(?i)\b((?:[a-z][\w-]+:(?:/{1,3}|[a-z0-9%])|www\d{0,3}[.]|[a-z0-9.\-]+[.][a-z]{2,4}/)(?:[^\s()<>\[\]]+|\(([^\s()<>\[\]]+|(\([^\s()<>\[\]]+\)))*\))+(?:\(([^\s()<>\[\]]+|(\([^\s()<>\[\]]+\)))*\)|[^\s`!(){};:'".,<>?\[\]]))"""  # noqa
)


def _is_url(url):
    return "://" in url and (url.startswith("http") or url.startswith("ftp"))


def _find_urls_in_text(text):
    return [url for url in HTTP_REGEX_PATTERN.findall(text) if _is_url(url[0])]


def construct_yaml_str(self, node):
    # Override the default string handling function
    # to always return unicode objects
    return self.construct_scalar(node)


yaml.Loader.add_constructor("tag:yaml.org,2002:str", construct_yaml_str)
yaml.SafeLoader.add_constructor("tag:yaml.org,2002:str", construct_yaml_str)


class ShedContext(NamedTuple):
    tsi: ToolShedInstance
    shed_config: dict
    config_owner: str

    def owner(self):
        owner = self.config_owner
        if owner is None:
            owner = username(self.tsi)
        return owner

    @property
    def label(self):
        return self.shed_config.get("label") or "tool shed"


def shed_init(ctx, path, **kwds):
    """Initialize a new shed repository."""
    if not os.path.exists(path):
        os.makedirs(path)
    shed_config_path = os.path.join(path, SHED_CONFIG_NAME)
    if not can_write_to_path(shed_config_path, **kwds):
        # .shed.yml exists and no --force sent.
        return 1

    create_failed = _create_shed_config(ctx, shed_config_path, **kwds)
    if create_failed:
        return 1

    repo_dependencies_path = os.path.join(path, REPO_DEPENDENCIES_CONFIG_NAME)
    from_workflow = kwds.get("from_workflow")

    if from_workflow:
        workflow_name = os.path.basename(from_workflow)
        workflow_target = os.path.join(path, workflow_name)
        if not os.path.exists(workflow_target):
            shutil.copyfile(from_workflow, workflow_target)

        if not can_write_to_path(repo_dependencies_path, **kwds):
            return 1

        repo_pairs = _parse_repos_from_workflow(from_workflow)
        repository_dependencies = RepositoryDependencies(repo_pairs)
        repository_dependencies.write_to_path(repo_dependencies_path)

    return 0


def install_arg_lists(ctx, paths, **kwds):
    """Build a list of install args for resolved repositories."""
    shed_context = get_shed_context(ctx, **kwds)
    install_args_list = []

    def process_repo(realized_repository):
        install_args_list.append(realized_repository.install_args(ctx, shed_context))
        return 0

    if not kwds.get("name"):
        exit_code = for_each_repository(ctx, process_repo, paths, **kwds)
        if exit_code:
            raise RuntimeError(PROBLEM_PROCESSING_REPOSITORY_MESSAGE)
    else:
        # Can only provide a single tool shed artifact to test
        name = kwds["name"]
        owner = kwds["owner"]
        changeset_revision = shed_context.tsi.repositories.get_ordered_installable_revisions(owner=owner, name=name)[-1]
        install_args_list.append(
            {
                "name": name,
                "owner": owner,
                "tool_shed_url": shed_context.tsi.base_url,
                "changeset_revision": changeset_revision,
            }
        )

    return install_args_list


def find_urls_for_xml(root):
    """Returns two lists: explicit package URLs, and help text URLs.

    For validating the user-facing URLs is it sensible to mimic
    a web browser user agent.
    """
    urls = []
    for packages in root.findall("package"):
        install_els = packages.findall("install")
        assert len(install_els) in (0, 1)

        if len(install_els) == 0:
            continue

        install_el = install_els[0]
        package = BasePackage(None, packages, install_el, readme=None)
        for action in package.get_all_actions():
            urls.extend([dl.url for dl in action.downloads()])

            for subaction in action.actions:
                if hasattr(subaction, "packages"):
                    urls.extend(subaction.packages)

    docs = []
    for help_text in root.findall("help"):
        for url in _find_urls_in_text(help_text.text):
            docs.append(url[0])

    return urls, docs


def handle_force_create(realized_repository, ctx, shed_context, **kwds):
    repo_id = realized_repository.find_repository_id(ctx, shed_context)
    if repo_id is None and kwds.get("force_repository_creation"):
        repo_id = realized_repository.create(ctx, shed_context)
    # failing to create the repo, give up
    return repo_id


def report_non_existent_repository(realized_repository):
    name = realized_repository.name
    error("Repository [%s] does not exist in the targeted Tool Shed." % name)
    return 2


def upload_repository(ctx, realized_repository, **kwds):
    """Upload a tool directory as a tarball to a tool shed."""
    path = realized_repository.path
    tar_path = kwds.get("tar")
    if not tar_path:
        tar_path = build_tarball(path, **kwds)
    if kwds.get("tar_only", False):
        name = realized_repository.pattern_to_file_name("shed_upload.tar.gz")
        shutil.copy(tar_path, name)
        return 0
    shed_context = get_shed_context(ctx, **kwds)
    update_kwds = {}
    _update_commit_message(ctx, realized_repository, update_kwds, **kwds)

    repo_id = handle_force_create(realized_repository, ctx, shed_context, **kwds)
    # failing to create the repo, give up
    if repo_id is None:
        return report_non_existent_repository(realized_repository)

    if kwds.get("check_diff", False):
        is_diff = diff_repo(ctx, realized_repository, **kwds) != 0
        if not is_diff:
            name = realized_repository.name
            info("Repository [%s] not different, skipping upload." % name)
            return 0

    # TODO: support updating repo information if it changes in the config file
    try:
        shed_context.tsi.repositories.update_repository(str(repo_id), tar_path, **update_kwds)
    except Exception as e:
        if isinstance(e, bioblend.ConnectionError) and e.status_code == 400 and '"No changes to repository."' in e.body:
            warn("Repository %s was not updated because there were no changes" % realized_repository.name)
            return 0
        message = api_exception_to_message(e)
        error("Could not update %s" % realized_repository.name)
        error(message)
        return -1
    info("Repository %s updated successfully." % realized_repository.name)
    return 0


def _update_commit_message(ctx, realized_repository, update_kwds, **kwds):
    message = kwds.get("message")
    git_rev = realized_repository.git_rev(ctx)
    git_repo = realized_repository.git_repo(ctx)
    if message is None:
        message = "planemo upload"
        if git_repo:
            message += " for repository %s" % git_repo
        if git_rev:
            message += " commit %s" % git_rev
    update_kwds["commit_message"] = message


def diff_repo(ctx, realized_repository, **kwds):
    """Compare two repositories (local or remote) and check for differences.

    Returns 0 if and only the repositories are effectively the same
    given supplied kwds for comparison description.
    """
    with temp_directory("tool_shed_diff_") as working:
        return _diff_in(ctx, working, realized_repository, **kwds)


def _diff_in(ctx, working, realized_repository, **kwds):
    path = realized_repository.path
    shed_target_source = kwds.get("shed_target_source")

    label_a = "_%s_" % (shed_target_source if shed_target_source else "workingdir")
    shed_target = kwds.get("shed_target", "B")
    if "/" in shed_target:
        shed_target = "custom_shed"
    label_b = "_%s_" % shed_target

    mine = os.path.join(working, label_a)
    other = os.path.join(working, label_b)

    shed_context = get_shed_context(ctx, read_only=True, **kwds)
    # In order to download the tarball, require repository ID...
    repo_id = realized_repository.find_repository_id(ctx, shed_context)
    if repo_id is None:
        error("shed_diff: Repository [%s] does not exist in the targeted Tool Shed." % realized_repository.name)
        # $ diff README.rst not_a_file 2&>1 /dev/null; echo $?
        # 2
        return 2
    info("Diffing repository [%s]" % realized_repository.name)
    download_tarball(
        ctx, shed_context, realized_repository, destination=other, clean=True, destination_is_pattern=False, **kwds
    )
    if shed_target_source:
        new_kwds = kwds.copy()
        new_kwds["shed_target"] = shed_target_source
        shed_context = get_shed_context(ctx, read_only=True, **new_kwds)
        download_tarball(
            ctx,
            shed_context,
            realized_repository,
            destination=mine,
            clean=True,
            destination_is_pattern=False,
            **new_kwds,
        )
    else:
        tar_path = build_tarball(path)
        os.mkdir(mine)
        shell(["tar", "-xzf", tar_path, "-C", mine])
        shutil.rmtree(tar_path, ignore_errors=True)

    output = kwds.get("output")
    raw = kwds.get("raw", False)
    xml_diff = 0
    if not raw:
        if output:
            with open(output, "w") as f:
                xml_diff = diff_and_remove(working, label_a, label_b, f)
        else:
            xml_diff = diff_and_remove(working, label_a, label_b, sys.stdout)

    cmd = ["diff", "-r", label_a, label_b]
    if output:
        with open(output, "ab") as fh:
            raw_diff = shell(cmd, cwd=working, stdout=fh)
    else:
        raw_diff = shell(cmd, cwd=working)
    exit = raw_diff or xml_diff
    if not raw:
        if xml_diff:
            ctx.vlog("One or more shed XML file(s) different!")
        if raw_diff:
            ctx.vlog("One or more non-shed XML file(s) different.")
        if not xml_diff and not raw_diff:
            ctx.vlog("No differences.")
    return exit


def shed_repo_config(ctx, path, name=None):
    shed_yaml_path = os.path.join(path, SHED_CONFIG_NAME)
    config = {}
    if os.path.exists(shed_yaml_path):
        with open(shed_yaml_path) as f:
            config = yaml.safe_load(f)

    if config is None:  # yaml may yield None
        config = {}
    _expand_raw_config(ctx, config, path, name=name)
    return config


def tool_shed_client(ctx=None, **kwds):
    return get_shed_context(ctx, **kwds).tsi


def get_shed_context(ctx=None, **kwds):
    read_only = kwds.get("read_only", False)
    shed_config, username = _shed_config_and_username(ctx, **kwds)

    def prop(key):
        return kwds.get("shed_%s" % key) or shed_config.get(key)

    url = _shed_config_to_url(shed_config)
    if read_only:
        key = None
        email = None
        password = None
    else:
        key = _find_shed_key(kwds, shed_config)
        email = prop("email")
        password = prop("password")

    tsi = tool_shed_instance(url, key, email, password)
    owner = username
    return ShedContext(tsi, shed_config, owner)


def tool_shed_url(ctx, **kwds):
    shed_config, _ = _shed_config_and_username(ctx, **kwds)
    return _shed_config_to_url(shed_config)


def _shed_config_and_username(ctx, **kwds):
    shed_target = kwds.get("shed_target")
    global_config = getattr(ctx, "global_config", {})
    if global_config and "sheds" in global_config:
        sheds_config = global_config["sheds"]
        shed_config = sheds_config.get(shed_target, {}) or {}
    else:
        shed_config = {}

    if "url" not in shed_config:
        if shed_target and shed_target in SHED_SHORT_NAMES:
            shed_config["url"] = SHED_SHORT_NAMES[shed_target]
        else:
            shed_config["url"] = shed_target

    if "label" not in shed_config:
        if shed_target and shed_target in SHED_LABELS:
            shed_config["label"] = SHED_LABELS[shed_target]
        else:
            shed_config["label"] = "custom tool shed at %s" % shed_target

    default_shed_username = global_config.get("shed_username")
    username = shed_config.get("username", default_shed_username)

    return shed_config, username


def _find_shed_key(kwds, shed_config):
    shed_key = kwds.get("shed_key")
    if shed_key is None:
        shed_key_from_env = kwds.get("shed_key_from_env")
        if shed_key_from_env is not None:
            shed_key = os.environ[shed_key_from_env]
    if shed_key is None:
        shed_key = shed_config.get("key")
    return shed_key


def find_repository_id(ctx, shed_context, path, **kwds):
    repo_config = kwds.get("config")
    if repo_config is None:
        name = kwds.get("name")
        repo_config = shed_repo_config(ctx, path, name=name)
    name = repo_config["name"]
    find_kwds = kwds.copy()
    if "name" in find_kwds:
        del find_kwds["name"]
    return _find_repository_id(ctx, shed_context, name, repo_config, **find_kwds)


def _find_repository_id(ctx, shed_context, name, repo_config, **kwds):
    # TODO: modify to consume shed_context
    owner = _owner(ctx, repo_config, shed_context, **kwds)
    matching_repository = find_repository(shed_context.tsi, owner, name)
    if matching_repository is None:
        raise Exception(f"Failed to find repository for owner/name {owner}/{name}")
    else:
        repo_id = matching_repository["id"]
        return repo_id


def _owner(ctx, repo_config, shed_context=None, **kwds):
    owner = kwds.get("owner") or repo_config.get("owner")
    if owner is None:
        if shed_context is None and "shed_target" in kwds:
            shed_context = get_shed_context(ctx, **kwds)
        if shed_context is not None:
            owner = shed_context.owner()
    return owner


def _expand_raw_config(ctx, config, path, name=None):
    name_input = name
    if "name" not in config:
        config["name"] = name
    if config["name"] is None:
        config["name"] = path_to_repo_name(path)

    default_include = config.get("include", ["**"])
    repos = config.get("repositories")
    auto_tool_repos = config.get("auto_tool_repositories", False)
    suite_config = config.get("suite", False)

    if repos and auto_tool_repos:
        raise Exception(AUTO_REPO_CONFLICT_MESSAGE)
    if auto_tool_repos and name_input:
        raise Exception(AUTO_NAME_CONFLICT_MESSAGE)
    if auto_tool_repos:
        repos = _build_auto_tool_repos(ctx, path, config, auto_tool_repos)
    if suite_config:
        if repos is None:
            repos = odict.odict()
        _build_suite_repo(config, repos, suite_config)
    # If repositories aren't defined, just define a single
    # one based on calculated name and including everything
    # by default.
    if repos is None:
        repos = {config["name"]: {"include": default_include}}
    config["repositories"] = repos


def _build_auto_tool_repos(ctx, path, config, auto_tool_repos):
    default_include = config.get("include", ["**"])
    tool_source_pairs = list(yield_tool_sources(ctx, path, recursive=True))
    paths = [_[0] for _ in tool_source_pairs]
    excludes = _shed_config_excludes(config)

    def _build_repository(tool_path, tool_source):
        tool_id = tool_source.parse_id().lower()
        tool_name = tool_source.parse_name()
        description = tool_source.parse_description()
        template_vars = dict(
            tool_id=tool_id,
            tool_name=tool_name,
            description=description,
        )
        other_paths = paths[:]
        other_paths.remove(tool_path)
        tool_excludes = excludes + list(other_paths)
        repo_dict = {
            "include": default_include,
            "exclude": tool_excludes,
        }
        for key in ["name", "description", "long_description"]:
            template_key = "%s_template" % key
            template = auto_tool_repos.get(template_key)
            if template:
                value = templates.render(template, **template_vars)
                repo_dict[key] = value
        return repo_dict

    repos = odict.odict()
    for tool_path, tool_source in tool_source_pairs:
        repository_config = _build_repository(tool_path, tool_source)
        repository_name = repository_config["name"]
        repos[repository_name] = repository_config
    return repos


def _build_suite_repo(config, repos, suite_config):
    name = suite_config.get("name")
    if not name:
        raise Exception("suite requires a 'name'.")
    description = suite_config.get("description", "")
    long_description = suite_config.get("long_description")
    owner = config["owner"]
    repo_type = suite_config.get("type", REPO_TYPE_SUITE)

    repo_pairs = [(repo_dict.get("owner") or owner, repo_name) for repo_name, repo_dict in repos.items()]
    extra_repos = suite_config.get("include_repositories", {})
    repo_pairs += [(_["owner"], _["name"]) for _ in extra_repos]

    repository_dependencies = RepositoryDependencies(repo_pairs, description)

    repo = {
        "_files": {REPO_DEPENDENCIES_CONFIG_NAME: str(repository_dependencies)},
        "include": [],
        "name": name,
        "description": description,
        "type": repo_type,
    }
    if long_description:
        repo["long_description"] = long_description
    repos[name] = repo


def update_repository_for(ctx, tsi, id, repo_config):
    name = repo_config["name"]
    description = repo_config.get("description")
    long_description = repo_config.get("long_description")
    remote_repository_url = repo_config.get("remote_repository_url")
    homepage_url = repo_config.get("homepage_url")
    categories = repo_config.get("categories", [])
    category_ids = find_category_ids(tsi, categories)

    _ensure_shed_description(description)

    repo = tsi.repositories.update_repository_metadata(
        id,
        name=name,
        synopsis=description,
        description=long_description,
        remote_repository_url=remote_repository_url,
        homepage_url=homepage_url,
        category_ids=category_ids,
    )
    return repo


def create_repository_for(ctx, tsi, name, repo_config):
    description = repo_config.get("description")
    long_description = repo_config.get("long_description")
    repo_type = shed_repo_type(repo_config, name)
    remote_repository_url = repo_config.get("remote_repository_url")
    homepage_url = repo_config.get("homepage_url")
    categories = repo_config.get("categories", [])
    category_ids = find_category_ids(tsi, categories)

    _ensure_shed_description(description)

    repo = tsi.repositories.create_repository(
        name=name,
        synopsis=description,
        description=long_description,
        type=repo_type,
        remote_repository_url=remote_repository_url,
        homepage_url=homepage_url,
        category_ids=category_ids,
    )
    return repo


def download_tarball(ctx, shed_context, realized_repository, **kwds):
    repo_id = realized_repository.find_repository_id(ctx, shed_context)
    if repo_id is None:
        message = "Unable to find repository id, cannot download."
        error(message)
        raise Exception(message)
    destination_pattern = kwds.get("destination", "shed_download.tar.gz")
    if kwds.get("destination_is_pattern", True):
        destination = realized_repository.pattern_to_file_name(destination_pattern)
    else:
        destination = destination_pattern
    to_directory = not destination.endswith("gz")
    download_tar(shed_context.tsi, repo_id, destination, to_directory=to_directory)
    if to_directory:
        clean = kwds.get("clean", False)
        if clean:
            archival_file = os.path.join(destination, ".hg_archival.txt")
            if os.path.exists(archival_file):
                os.remove(archival_file)


def build_tarball(realized_path, **kwds):
    """Build a tool-shed tar ball for the specified path, caller is
    responsible for deleting this file.
    """

    # Simplest solution to sorting the files is to use a list,
    files = []
    for dirpath, _dirnames, filenames in os.walk(realized_path):
        for f in filenames:
            files.append(os.path.join(dirpath, f))
    files.sort()

    fd, temp_path = mkstemp()
    try:
        tar = tarfile.open(temp_path, "w:gz", dereference=True)
        try:
            for raw in files:
                name = os.path.relpath(raw, realized_path)
                tar.add(os.path.join(realized_path, name), arcname=name)
        finally:
            tar.close()
    finally:
        os.close(fd)
    return temp_path


def find_raw_repositories(ctx, paths, **kwds):
    """Return a list of "raw" repository objects for each repo on paths."""
    raw_repo_objects = []
    for path in paths:
        raw_repo_objects.extend(_find_raw_repositories(ctx, path, **kwds))
    return raw_repo_objects


def for_each_repository(ctx, function, paths, **kwds):
    ret_codes = []
    for path in paths:
        with _path_on_disk(ctx, path) as raw_path:
            try:
                for realized_repository in _realize_effective_repositories(ctx, raw_path, **kwds):
                    ret_codes.append(function(realized_repository))
            except RealizationException:
                error(REALIZAION_PROBLEMS_MESSAGE)
                return 254

    return coalesce_return_codes(ret_codes)


def path_to_repo_name(path):
    return os.path.basename(os.path.abspath(path))


def shed_repo_type(config, name):
    repo_type = config.get("type")
    if repo_type is None:
        if name.startswith("package_"):
            repo_type = REPO_TYPE_TOOL_DEP
        elif name.startswith("suite_"):
            repo_type = REPO_TYPE_SUITE
        else:
            repo_type = REPO_TYPE_UNRESTRICTED
    return repo_type


def _shed_config_to_url(shed_config):
    url = shed_config["url"]
    if not url.startswith("http"):
        message = f"Invalid shed url specified [{url}]. Please specify a valid HTTP address or one of {list(SHED_SHORT_NAMES.keys())}"
        raise ValueError(message)
    return url


def _realize_effective_repositories(ctx, path, **kwds):
    """Expands folders in a source code repository into tool shed
    repositories.

    Each folder may have nested repositories and each folder may corresponding
    to many repositories (for instance if a folder has n tools in the source
    code repository but are published to the tool shed as one repository per
    tool).
    """
    raw_repo_objects = _find_raw_repositories(ctx, path, **kwds)
    failed = False
    with temp_directory() as base_dir:
        for raw_repo_object in raw_repo_objects:
            if isinstance(raw_repo_object, Exception):
                _handle_realization_error(raw_repo_object, **kwds)
                failed = True
                continue

            realized_repos = raw_repo_object.realizations(ctx, base_dir, **kwds)
            for realized_repo in realized_repos:
                if isinstance(realized_repo, Exception):
                    _handle_realization_error(realized_repo, **kwds)
                    failed = True
                    continue
                yield realized_repo
    if failed:
        raise RealizationException()


def _create_shed_config(ctx, path, **kwds):
    name = kwds.get("name") or path_to_repo_name(os.path.dirname(path))
    name_invalid = validate_repo_name(name)
    if name_invalid:
        error(name_invalid)
        return 1

    owner = kwds.get("owner")
    if owner is None:
        owner = ctx.global_config.get("shed_username")
    owner_invalid = validate_repo_owner(owner)
    if owner_invalid:
        error(owner_invalid)
        return 1
    description = kwds.get("description") or name
    long_description = kwds.get("long_description")
    remote_repository_url = kwds.get("remote_repository_url")
    homepage_url = kwds.get("homepage_url")
    categories = kwds.get("category", [])
    config = dict(
        name=name,
        owner=owner,
        description=description,
        long_description=long_description,
        remote_repository_url=remote_repository_url,
        homepage_url=homepage_url,
        categories=categories,
    )
    # Remove empty entries...
    for k in list(config.keys()):
        if config[k] is None:
            del config[k]

    with open(path, "w") as f:
        yaml.safe_dump(config, f)


def _parse_repos_from_workflow(path):
    with open(path) as f:
        workflow_json = json.load(f)
    steps = workflow_json["steps"]
    tool_ids = set()
    for value in steps.values():
        step_type = value["type"]
        if step_type != "tool":
            continue
        tool_id = value["tool_id"]
        if "/repos/" in tool_id:
            tool_ids.add(tool_id)

    repo_pairs = set()
    for tool_id in tool_ids:
        tool_repo_info = tool_id.split("/repos/", 1)[1]
        tool_repo_parts = tool_repo_info.split("/")
        owner = tool_repo_parts[0]
        name = tool_repo_parts[1]
        repo_pairs.add((owner, name))

    return repo_pairs


@contextlib.contextmanager
def _path_on_disk(ctx, path):
    git_path = None
    if path.startswith("git:"):
        git_path = path
    elif path.startswith("git+"):
        git_path = path[len("git+") :]
    if git_path is None:
        yield path
    else:
        with temp_directory() as git_repo:
            git.clone(ctx, git_path, git_repo)
            yield git_repo


def _find_raw_repositories(ctx, path, **kwds):
    name = kwds.get("name")
    recursive = kwds.get("recursive", False)

    shed_file_dirs = []
    for pattern in REPO_METADATA_FILES:
        shed_file_dirs.extend(find_matching_directories(path, pattern, recursive=recursive))
    config_name = None
    if len(shed_file_dirs) == 1:
        shed_file_dir = shed_file_dirs[0]
        try:
            config = shed_repo_config(ctx, shed_file_dir, name=name)
        except Exception as e:
            error_message = PARSING_PROBLEM % (shed_file_dir, e)
            exception = RuntimeError(error_message)
            _handle_realization_error(exception, **kwds)
            return [exception]
        config_name = config.get("name")

    if len(shed_file_dirs) > 1 and name is not None:
        raise Exception(NAME_INVALID_MESSAGE)
    if config_name is not None and name is not None:
        if config_name != name:
            raise Exception(CONFLICTING_NAMES_MESSAGE)
    raw_dirs = shed_file_dirs or [path]
    kwds_copy = kwds.copy()
    kwds_copy["name"] = name
    return _build_raw_repo_objects(ctx, raw_dirs, **kwds_copy)


def _build_raw_repo_objects(ctx, raw_dirs, **kwds):
    """
    From specific directories with .shed.yml files or specified directly from
    the command-line build abstract description of directories that should be
    expanded out into shed repositories.
    """
    multiple = len(raw_dirs) > 1
    name = kwds.get("name")

    # List of RawRepositoryDirectories or parsing failures if
    # fail_fast is not enabled.
    raw_repo_objects = []
    for raw_dir in raw_dirs:
        try:
            config = shed_repo_config(ctx, raw_dir, name=name)
        except Exception as e:
            error_message = PARSING_PROBLEM % (raw_dir, e)
            exception = RuntimeError(error_message)
            _handle_realization_error(exception, **kwds)
            raw_repo_objects.append(exception)
            continue
        raw_repo_object = RawRepositoryDirectory(raw_dir, config, multiple)
        raw_repo_objects.append(raw_repo_object)
    return raw_repo_objects


class RepositoryDependencies:
    """Abstraction for shed repository_dependencies.xml files."""

    def __init__(self, repo_pairs, description=None):
        self.repo_pairs = repo_pairs
        self.description = description or ""

    def __str__(self):
        contents = '<repositories description="%s">' % self.description
        line_template = '  <repository owner="%s" name="%s" />\n'
        for owner, name in self.repo_pairs:
            contents += line_template % (owner, name)
        contents += "</repositories>"
        return contents

    def write_to_path(self, path):
        with open(path, "w") as f:
            f.write(str(self))


class RawRepositoryDirectory:
    def __init__(self, path, config, multiple):
        self.path = path
        self.config = config
        self.name = config["name"]
        self.type = shed_repo_type(config, self.name)
        self.multiple = multiple  # operation over many repos?

    def _hash(self, name):
        return hashlib.md5(name.encode("utf-8")).hexdigest()

    def realizations(self, ctx, parent_directory, **kwds):
        names = self._repo_names()

        for name in names:
            directory = os.path.join(parent_directory, self._hash(name), name)
            multiple = self.multiple or len(names) > 1
            if not os.path.exists(directory):
                os.makedirs(directory)
            r_kwds = kwds.copy()
            if "name" in r_kwds:
                del r_kwds["name"]
            yield self._realize_to(ctx, directory, name, multiple, **r_kwds)

    def _realize_to(self, ctx, directory, name, multiple, **kwds):
        fail_on_missing = kwds.get("fail_on_missing", True)
        ignore_list = []
        config = self._realize_config(name)
        config["owner"] = _owner(ctx, config, **kwds)

        excludes = _shed_config_excludes(config)
        for exclude in excludes:
            ignore_list.extend(_glob(self.path, exclude))

        realized_files = self._realized_files(name)
        missing = realized_files.include_failures
        if missing and fail_on_missing:
            msg = "Failed to include files for %s" % missing
            return RuntimeError(msg)

        for realized_file in realized_files.files:
            relative_dest = realized_file.dest
            implicit_ignore = self._implicit_ignores(relative_dest)
            explicit_ignore = realized_file.absolute_src in ignore_list
            if implicit_ignore or explicit_ignore:
                continue
            realized_file.realize_to(directory)

        for name, contents in config.get("_files", {}).items():
            path = os.path.join(directory, name)
            with open(path, "w") as f:
                f.write(contents)

        return RealizedRepository(
            realized_path=directory,
            real_path=self.path,
            config=config,
            multiple=multiple,
            missing=missing,
        )

    def _repo_names(self):
        return self.config.get("repositories").keys()

    def _realized_files(self, name):
        config = self._realize_config(name)
        realized_files = []
        missing = []
        for include_info in config["include"]:
            if not isinstance(include_info, dict):
                include_info = {"source": include_info}
            source_list = include_info.get("source")
            if not isinstance(source_list, list):
                source_list = [source_list]
            # Preprocess any entries with a source list into copies
            # with a single source entry:
            for source in source_list:
                include = include_info.copy()
                include["source"] = source
                included = RealizedFile.realized_files_for(self.path, include)
                if not included:
                    missing.append(include)
                else:
                    realized_files.extend(included)
        return RealizedFiles(realized_files, missing)

    def _realize_config(self, name):
        config = copy.deepcopy(self.config)
        config["name"] = name
        repo_config = config.get("repositories", {}).get(name, {})
        config.update(repo_config)
        if "repositories" in config:
            del config["repositories"]
        return config

    def _implicit_ignores(self, relative_path):
        # Filter out "unwanted files" :) like READMEs for special
        # repository types.
        if self.type == REPO_TYPE_TOOL_DEP:
            if relative_path != TOOL_DEPENDENCIES_CONFIG_NAME:
                return True

        if self.type == REPO_TYPE_SUITE:
            if relative_path != REPO_DEPENDENCIES_CONFIG_NAME:
                return True

        name = os.path.basename(relative_path)
        for dvcs_prefix in [".git", ".hg"]:
            if relative_path.startswith(dvcs_prefix):
                return True

        if name.startswith(".svn"):
            return True

        for pattern in PLANEMO_FILES:
            if fnmatch.fnmatch(name, pattern):
                return True
        return False


class RealizedFiles(NamedTuple):
    files: list
    include_failures: list


class RealizedFile:
    def __init__(self, src_root, src, dest):
        """Create object mapping from file system to tar-ball.

        * src_root - source root (i.e. folder with .shed.yml file)
        * src - location of source file, relative to src_root
        * dest - destination path, relative to root of tar-ball.
        """
        if dest == ".":
            raise ValueError("Destination for %r should be a full filename!" % src)
        self.src_root = src_root
        self.src = src
        self.dest = dest

    @property
    def absolute_src(self):
        return os.path.abspath(os.path.join(self.src_root, self.src))

    def realize_to(self, directory):
        source_path = self.absolute_src
        if os.path.islink(source_path):
            source_path = os.path.realpath(source_path)
        relative_dest = self.dest
        assert relative_dest != "."
        target_path = os.path.join(directory, relative_dest)
        target_exists = os.path.exists(target_path)
        # info("realize_to %r --> %r" % (source_path, target_path))
        if not target_exists:
            target_dir = os.path.dirname(target_path)
            if not os.path.exists(target_dir):
                os.makedirs(target_dir)
            if os.path.isdir(source_path):
                os.makedirs(target_path)
            else:
                os.symlink(source_path, target_path)

    @staticmethod
    def realized_files_for(path, include_info):
        if not isinstance(include_info, dict):
            include_info = {"source": include_info}
        source = include_info.get("source")
        abs_source = os.path.join(path, source)
        destination = include_info.get("destination")
        strip_components = include_info.get("strip_components", 0)
        if destination is None:
            destination = "./"
        if not destination.endswith("/"):
            # Check if source using wildcards (directory gets implicit wildcard)
            # Should we use a regular exoression to catch [A-Z] style patterns?
            if "*" in source or "?" in source or os.path.isdir(abs_source):
                raise ValueError(
                    "destination must be a directory (with trailing slash) if source is a folder or uses wildcards"
                )
        realized_files = []
        for globbed_file in _glob(path, source):
            src = os.path.relpath(globbed_file, path)
            if not destination.endswith("/"):
                # Given a filename, just use it!
                dest = destination
                if strip_components:
                    raise ValueError("strip_components should not be used if destination is a filename")
            else:
                # Destination is a directory...
                if not strip_components:
                    dest = src
                elif "/../" in globbed_file:
                    # Can't work from src=os.path.relpath(globbed_file, path) as lost any '..'
                    assert globbed_file.startswith(path + "/")
                    dest = "/".join(globbed_file[len(path) + 1 :].split("/")[strip_components:])
                else:
                    dest = "/".join(src.split("/")[strip_components:])
                # Now apply the specified output directory:
                dest = os.path.join(destination, dest)
            realized_files.append(RealizedFile(path, src, os.path.normpath(dest)))
        return realized_files

    def __str__(self):
        return f"RealizedFile[src={self.src},dest={self.dest},src_root={self.src_root}]"


class RealizedRepository:
    def __init__(self, realized_path, real_path, config, multiple, missing):
        self.path = realized_path
        self.real_path = real_path
        self.config = config
        self.name = config["name"]
        self.multiple = multiple
        self.missing = missing

    @property
    def owner(self):
        return self.config["owner"]

    @property
    def repository_type(self):
        return shed_repo_type(self.config, self.name)

    @property
    def is_package(self):
        return self.repository_type == REPO_TYPE_TOOL_DEP

    @property
    def is_suite(self):
        return self.repository_type == REPO_TYPE_SUITE

    @property
    def repo_dependencies_path(self):
        return os.path.join(self.path, REPO_DEPENDENCIES_CONFIG_NAME)

    @property
    def tool_dependencies_path(self):
        return os.path.join(self.path, TOOL_DEPENDENCIES_CONFIG_NAME)

    def git_rev(self, ctx):
        return git.rev_if_git(ctx, self.real_path)

    def git_repo(self, ctx):
        return self.config.get("remote_repository_url")

    def pattern_to_file_name(self, pattern):
        if not self.multiple:
            return pattern

        name = self.config["name"]
        suffix = "_%s" % name.replace("-", "_")

        if "." not in pattern:
            return pattern + suffix
        else:
            parts = pattern.split(".", 1)
            return parts[0] + suffix + "." + parts[1]

    def find_repository_id(self, ctx, shed_context):
        try:
            repo_id = _find_repository_id(
                ctx,
                shed_context,
                name=self.name,
                repo_config=self.config,
            )
            return repo_id
        except Exception as e:
            message = api_exception_to_message(e)
            error("Could not update %s" % self.name)
            error(message)
        return None

    def create(self, ctx, shed_context):
        """Wrapper for creating the endpoint if it doesn't exist"""
        context_owner = shed_context.owner()
        config_owner = self.config.get("owner")
        if context_owner and config_owner and context_owner != config_owner:
            # This is broken because context_owner is incorrect if using an API key.
            # message = INCORRECT_OWNER_MESSAGE % (config_owner, context_owner)
            # raise Exception(message)
            pass

        def _create():
            repo = create_repository_for(
                ctx,
                shed_context.tsi,
                self.name,
                self.config,
            )
            return repo["id"]

        return self._with_ts_exception_handling(_create)

    def update(self, ctx, shed_context, id):
        """Wrapper for update the repository metadata."""

        def _update():
            repo = update_repository_for(
                ctx,
                shed_context.tsi,
                id,
                self.config,
            )
            return repo

        return self._with_ts_exception_handling(_update)

    def _with_ts_exception_handling(self, f):
        try:
            return f()
        except Exception as e:
            # TODO: galaxyproject/bioblend#126
            try:
                upstream_error = json.loads(e.read())
                error(upstream_error["err_msg"])
            except Exception:
                error(unicodify(e))
            return None

    def latest_installable_revision(self, ctx: "PlanemoCliContext", shed_context: ShedContext):
        repository_id = self.find_repository_id(ctx, shed_context)
        return latest_installable_revision(shed_context.tsi, repository_id)

    def install_args(self, ctx, shed_context):
        """Arguments for bioblend's install_repository_revision
        to install this repository against supplied tsi.
        """
        tool_shed_url = shed_context.tsi.base_url
        return dict(
            tool_shed_url=tool_shed_url,
            name=self.name,
            owner=self.owner,
            changeset_revision=self.latest_installable_revision(ctx, shed_context),
        )


def _glob(path, pattern):
    pattern = os.path.join(path, pattern)
    if os.path.isdir(pattern):
        pattern = "%s/**" % pattern
    return glob.glob(pattern)


def _shed_config_excludes(config):
    return config.get("ignore", []) + config.get("exclude", [])


def _handle_realization_error(exception, **kwds):
    fail_fast = kwds.get("fail_fast", False)
    if fail_fast:
        raise exception
    else:
        error(unicodify(exception))


def _ensure_shed_description(description):
    # description is required, as is name.
    if description is None:
        message = "description required for automatic creation or update of shed metadata."
        raise ValueError(message)


def validate_repo_name(name):
    def _build_error(descript):
        return f"Repository name [{name}] invalid. {descript}"

    msg = None
    if len(name) < 2:
        msg = _build_error("Repository names must be at least 2 characters in length.")
    if len(name) > 80:
        msg = _build_error("Repository names cannot be more than 80 characters in length.")
    if not VALID_REPOSITORYNAME_RE.match(name):
        msg = _build_error("Repository names must contain only lower-case letters, numbers and underscore.")
    return msg


def validate_repo_owner(owner):
    def _build_error(descript):
        return f"Owner [{owner}] invalid. {descript}"

    msg = None
    if len(owner) < 3:
        msg = _build_error("Owner must be at least 3 characters in length")
    if len(owner) > 255:
        msg = _build_error("Owner cannot be more than 255 characters in length")
    if not (VALID_PUBLICNAME_RE.match(owner)):
        msg = _build_error("Owner must contain only lower-case letters, numbers, dots, underscores, and '-'")
    return msg


class RealizationException(Exception):
    """This exception indicates there was a problem while
    realizing effective repositories for a shed command. As a
    precondition - the user has already been informed with error().
    """


__all__ = (
    "api_exception_to_message",
    "CURRENT_CATEGORIES",
    "diff_repo",
    "download_tarball",
    "find_raw_repositories",
    "for_each_repository",
    "get_shed_context",
    "path_to_repo_name",
    "REPO_TYPE_SUITE",
    "REPO_TYPE_TOOL_DEP",
    "REPO_TYPE_UNRESTRICTED",
    "shed_init",
    "tool_shed_client",  # Deprecated...
    "tool_shed_url",
)
