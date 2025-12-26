"""Module contianing the :class:`GalaxyEngine` implementation of :class:`Engine`."""

import abc
import contextlib
from typing import (
    Callable,
    List,
    Optional,
    TYPE_CHECKING,
)

from galaxy.tool_util.verify import interactor

from planemo import io
from planemo.database.postgres_singularity import SingularityPostgresDatabaseSource
from planemo.galaxy.activity import (
    execute,
    execute_rerun,
    GalaxyBaseRunResponse,
)
from planemo.galaxy.config import external_galaxy_config
from planemo.galaxy.serve import serve_daemon
from planemo.runnable import (
    DelayedGalaxyToolTestCase,
    ExternalGalaxyToolTestCase,
    GALAXY_TOOLS_PREFIX,
    Rerunnable,
    RunnableType,
)
from .interface import BaseEngine

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

INSTALLING_MESSAGE = "Installing repositories - this may take some time..."


class GalaxyEngine(BaseEngine, metaclass=abc.ABCMeta):
    """An :class:`Engine` implementation backed by a managed Galaxy.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    handled_runnable_types = [
        RunnableType.cwl_tool,
        RunnableType.cwl_workflow,
        RunnableType.galaxy_workflow,
        RunnableType.galaxy_tool,
        RunnableType.galaxy_datamanager,
        RunnableType.directory,
    ]

    def _run(self, runnables, job_paths, output_collectors: Optional[List[Callable]] = None):
        """Run job in Galaxy."""
        results = []
        if not output_collectors:
            output_collectors = [lambda x: None] * len(runnables)

        with self.ensure_runnables_served(runnables) as config:
            if self._ctx.verbose:
                self._ctx.log(f"Running Galaxy with API configuration [{config.user_api_config}]")
            for runnable, job_path, collect_output in zip(runnables, job_paths, output_collectors):
                self._ctx.vlog(f"Serving artifact [{runnable}] with Galaxy.")
                self._ctx.vlog(f"Running job path [{job_path}]")
                run_response = execute(self._ctx, config, runnable, job_path, **self._kwds)
                results.append(run_response)
                if collect_output is not None:
                    collect_output(run_response)

        return results

    @abc.abstractmethod
    def ensure_runnables_served(self, runnables):
        """Use a context manager and describe Galaxy instance with runnables being served."""

    def _run_test_cases(self, test_cases, test_timeout):
        test_results = []
        file_based_test_cases = []
        embedded_test_cases = []
        # TODO: unify interface so we don't need to split test cases
        for test_case in test_cases:
            if isinstance(test_case, ExternalGalaxyToolTestCase):
                embedded_test_cases.append(test_case)
            else:
                file_based_test_cases.append(test_case)
        if file_based_test_cases:
            test_results.extend(super()._run_test_cases(file_based_test_cases, test_timeout))
        if embedded_test_cases:
            runnables = [test_case.runnable for test_case in embedded_test_cases]
            with self.ensure_runnables_served(runnables) as config:
                expanded_test_cases = expand_test_cases(config, embedded_test_cases)
                for test_case in expanded_test_cases:
                    galaxy_interactor_kwds = {
                        "galaxy_url": config.galaxy_url,
                        "master_api_key": config.master_api_key,
                        "api_key": config.user_api_key,
                        "keep_outputs_dir": self._kwds.get("test_data_target_dir"),
                    }
                    tool_id = test_case.tool_id
                    test_index = test_case.test_index
                    tool_version = test_case.tool_version
                    galaxy_interactor = interactor.GalaxyInteractorApi(**galaxy_interactor_kwds)

                    def _register_job_data(job_data):
                        test_results.append(
                            {
                                "id": tool_id + "-" + str(test_index),
                                "has_data": True,
                                "data": job_data,
                            }
                        )

                    verbose = self._ctx.verbose
                    try:
                        if verbose:
                            # TODO: this is pretty hacky, it'd be better to send a stream
                            # and capture the output information somehow.
                            interactor.VERBOSE_GALAXY_ERRORS = True

                        interactor.verify_tool(
                            tool_id,
                            galaxy_interactor,
                            test_index=test_index,
                            tool_version=tool_version,
                            register_job_data=_register_job_data,
                            maxseconds=test_timeout,
                            quiet=not verbose,
                        )
                    except Exception:
                        pass

        return test_results


class LocalManagedGalaxyEngine(GalaxyEngine):
    """An :class:`Engine` implementation backed by a managed Galaxy.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    @contextlib.contextmanager
    def ensure_runnables_served(self, runnables):
        # TODO: define an interface for this - not everything in config would make sense for a
        # pre-existing Galaxy interface.
        with serve_daemon(self._ctx, runnables, **self._serve_kwds()) as config:
            if "install_args_list" in self._serve_kwds():
                self.shed_install(config)
            yield config

    def shed_install(self, config):
        kwds = self._serve_kwds()
        install_args_list = kwds["install_args_list"]
        install_deps = not kwds.get("skip_dependencies", False)
        print(INSTALLING_MESSAGE)
        io.info(INSTALLING_MESSAGE)
        for install_args in install_args_list:
            install_args["install_tool_dependencies"] = install_deps
            install_args["install_repository_dependencies"] = True
            install_args["new_tool_panel_section_label"] = "Shed Installs"
            config.install_repo(**install_args)
        try:
            config.wait_for_all_installed()
        except Exception:
            if self._ctx.verbose:
                print("Failed to install tool repositories, Galaxy log:")
                print(config.log_contents)
                print("Galaxy root:")
                io.shell(["ls", config.galaxy_root])
            raise

    def _serve_kwds(self):
        return self._kwds.copy()


class LocalManagedGalaxyEngineWithSingularityDB(LocalManagedGalaxyEngine):
    def run(self, runnables, job_paths, output_collectors: Optional[List[Callable]] = None):
        with SingularityPostgresDatabaseSource(**self._kwds.copy()):
            run_responses = super().run(runnables, job_paths, output_collectors)
        return run_responses


class DockerizedManagedGalaxyEngine(LocalManagedGalaxyEngine):
    """An :class:`Engine` implementation backed by Galaxy running in Docker.

    More information on Galaxy can be found at http://galaxyproject.org/.
    """

    def _serve_kwds(self):
        serve_kwds = self._kwds.copy()
        serve_kwds["dockerize"] = True
        return serve_kwds


class ExternalGalaxyEngine(GalaxyEngine):
    """An :class:`Engine` implementation backed by an external Galaxy instance."""

    @contextlib.contextmanager
    def ensure_runnables_served(self, runnables):
        # TODO: ensure tools are available
        with external_galaxy_config(self._ctx, runnables, **self._kwds) as config:
            config.install_workflows()
            yield config

    def rerun(self, ctx: "PlanemoCliContext", rerunnable: Rerunnable, **kwds) -> GalaxyBaseRunResponse:
        with self.ensure_runnables_served([]) as config:
            rerun_response = execute_rerun(ctx, config, rerunnable, **kwds)
            return rerun_response


def expand_test_cases(config, test_cases):
    expanded_test_cases = []
    for test_case in test_cases:
        if not isinstance(test_case, DelayedGalaxyToolTestCase):
            expanded_test_cases.append(test_case)
        else:
            runnable = test_case.runnable
            tool_id = runnable.uri.split(GALAXY_TOOLS_PREFIX)[1]
            test_data = config.gi.tools._get(f"{tool_id}/test_data")
            for test_dict in test_data:
                expanded_test_cases.append(
                    ExternalGalaxyToolTestCase(
                        runnable,
                        tool_id=tool_id,
                        tool_version=test_dict["tool_version"],
                        test_index=test_dict["test_index"],
                        test_dict=test_dict,
                    )
                )
    return expanded_test_cases


__all__ = (
    "DockerizedManagedGalaxyEngine",
    "ExternalGalaxyEngine",
    "LocalManagedGalaxyEngine",
    "LocalManagedGalaxyEngineWithSingularityDB",
)
