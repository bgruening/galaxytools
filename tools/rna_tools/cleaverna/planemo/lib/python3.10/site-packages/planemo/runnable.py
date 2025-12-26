"""Describe artifacts that can be run, tested, and linted."""

import abc
import os
from enum import (
    auto,
    Enum,
)
from pathlib import Path
from typing import (
    Any,
    Dict,
    Iterable,
    List,
    NamedTuple,
    Optional,
    Union,
)
from urllib.parse import urlparse

import yaml
from bioblend.galaxy import GalaxyInstance
from galaxy.tool_util.cwl.parser import workflow_proxy
from galaxy.tool_util.loader_directory import (
    is_a_yaml_with_class,
    looks_like_a_cwl_artifact,
    looks_like_a_data_manager_xml,
    looks_like_a_tool_cwl,
    looks_like_a_tool_xml,
)
from galaxy.tool_util.parser import get_tool_source

from planemo.exit_codes import (
    EXIT_CODE_UNKNOWN_FILE_TYPE,
    ExitCodeException,
)
from planemo.galaxy.workflows import (
    describe_outputs,
    GALAXY_WORKFLOW_INSTANCE_PREFIX,
    GALAXY_WORKFLOWS_PREFIX,
    WorkflowOutput,
)
from planemo.io import error
from planemo.shed import DOCKSTORE_REGISTRY_CONF
from planemo.test import (
    check_output,
    for_collections,
)
from planemo.tools import yield_tool_sources_on_paths

TEST_SUFFIXES = ["-tests", "_tests", "-test", "_test"]
TEST_EXTENSIONS = [".yml", ".yaml", ".json"]

TEST_FILE_NOT_LIST_MESSAGE = "Invalid test definition file [%s] - file must contain a list of tests"
TEST_FIELD_MISSING_MESSAGE = "Invalid test definition [test #%d in %s] -defintion must field [%s]."
GALAXY_TOOLS_PREFIX = "gxid://tools/"


class RunnableType(Enum):
    galaxy_tool = auto()
    galaxy_datamanager = auto()
    galaxy_workflow = auto()
    cwl_tool = auto()
    cwl_workflow = auto()
    directory = auto()

    @property
    def has_tools(runnable_type):
        return runnable_type.name in ["galaxy_tool", "galaxy_datamanager", "cwl_tool", "directory"]

    @property
    def is_single_artifact(runnable_type):
        return runnable_type.name not in ["directory"]

    @property
    def test_data_in_parent_dir(runnable_type):
        return runnable_type.name in ["galaxy_datamanager"]

    @property
    def is_galaxy_artifact(runnable_type) -> bool:
        return "galaxy" in runnable_type.name

    @property
    def is_cwl_artifact(runnable_type) -> bool:
        return "cwl" in runnable_type.name


class Runnable(NamedTuple):
    """Abstraction describing tools and workflows."""

    uri: str
    type: RunnableType

    @property
    def path(self) -> str:
        uri = self.uri
        if self.is_remote_workflow_uri:
            parse_result = urlparse(uri)
            query = parse_result.query
            if query:
                assert query.startswith("runnable_path=")
                return query[len("runnable_path=") :]
            else:
                raise ValueError(f"Runnable with URI {uri} is remote resource without local path")
        else:
            return uri

    @property
    def has_path(self):
        try:
            self.path
            return True
        except ValueError:
            return False

    @property
    def is_remote_workflow_uri(self) -> bool:
        return self.uri.startswith((GALAXY_WORKFLOWS_PREFIX, GALAXY_WORKFLOW_INSTANCE_PREFIX))

    @property
    def test_data_search_path(self) -> str:
        """During testing, path to search for test data files."""
        if self.type.name in ["galaxy_datamanager"]:
            return os.path.join(os.path.dirname(self.path), os.path.pardir)
        else:
            return self.path

    @property
    def tool_data_search_path(self) -> str:
        """During testing, path to search for Galaxy tool data tables."""
        return self.test_data_search_path

    @property
    def data_manager_conf_path(self) -> Optional[str]:
        """Path of a Galaxy data manager configuration for runnable or None."""
        if self.type.name in ["galaxy_datamanager"]:
            return os.path.join(os.path.dirname(self.path), os.pardir, "data_manager_conf.xml")
        return None

    @property
    def has_tools(self) -> property:
        """Boolean indicating if this runnable corresponds to one or more tools."""
        return _runnable_delegate_attribute("has_tools")

    @property
    def is_single_artifact(self) -> property:
        """Boolean indicating if this runnable is a single artifact.

        Currently only directories are considered not a single artifact.
        """
        return _runnable_delegate_attribute("is_single_artifact")


class Rerunnable(NamedTuple):
    """Abstraction describing artifacts (histories, invocation, jobs) on external Galaxy instances with associated rerunnable and remappable jobs."""

    rerunnable_id: str
    rerunnable_type: str
    server_url: str


def _runnable_delegate_attribute(attribute: str) -> property:
    def getter(runnable):
        return getattr(runnable.type, attribute)

    return property(getter)


def workflows_from_dockstore_yaml(path):
    workflows = []
    parent_dir = Path(path).absolute().parent
    with open(path) as y:
        for workflow in yaml.safe_load(y).get("workflows", []):
            workflow_path = workflow.get("primaryDescriptorPath")
            if workflow_path:
                if workflow_path.startswith("/"):
                    workflow_path = workflow_path[1:]
            workflows.append(parent_dir.joinpath(workflow_path))
    return workflows


def workflow_dir_runnables(path: str) -> List[Runnable]:
    dockstore_path = os.path.join(path, DOCKSTORE_REGISTRY_CONF)
    if os.path.exists(dockstore_path):
        return [
            Runnable(str(path), RunnableType.galaxy_workflow) for path in workflows_from_dockstore_yaml(dockstore_path)
        ]
    return []


def tool_dir_runnables(path: str) -> List[Runnable]:
    return for_paths(tool_path for tool_path, _ in yield_tool_sources_on_paths(ctx=None, paths=[path]))


def for_path(path: str) -> Union[Runnable, List[Runnable]]:
    """Produce a class:`Runnable` for supplied path."""
    runnable_type = None
    if os.path.isdir(path):
        runnable = workflow_dir_runnables(path) or tool_dir_runnables(path)
        if runnable:
            return runnable
        runnable_type = RunnableType.directory
    elif looks_like_a_tool_cwl(path):
        runnable_type = RunnableType.cwl_tool
    elif looks_like_a_data_manager_xml(path):
        runnable_type = RunnableType.galaxy_datamanager
    elif looks_like_a_tool_xml(path):
        runnable_type = RunnableType.galaxy_tool
    elif is_a_yaml_with_class(path, ["GalaxyWorkflow"]):
        runnable_type = RunnableType.galaxy_workflow
    elif path.endswith(".ga"):
        runnable_type = RunnableType.galaxy_workflow
    elif looks_like_a_cwl_artifact(path, ["Workflow"]):
        runnable_type = RunnableType.cwl_workflow
    else:
        # Check to see if it is a Galaxy workflow with a different extension
        try:
            with open(path) as f:
                as_dict = yaml.safe_load(f)
            if as_dict.get("a_galaxy_workflow", False):
                runnable_type = RunnableType.galaxy_workflow
        except Exception:
            pass

    if runnable_type is None:
        error(f"Unable to determine runnable type for path [{path}]")
        raise ExitCodeException(EXIT_CODE_UNKNOWN_FILE_TYPE)

    return Runnable(path, runnable_type)


def for_paths(paths: Iterable[str]) -> List[Runnable]:
    """Return a specialized list of Runnable objects for paths."""
    runnables = []
    for path in paths:
        runnables_for_path = for_path(path)
        if isinstance(runnables_for_path, list):
            runnables.extend(runnables_for_path)
        else:
            runnables.append(runnables_for_path)
    return runnables


def for_uri(uri: str) -> Runnable:
    """Produce a class:`Runnable` for supplied Galaxy workflow or tool ID."""
    runnable_type = RunnableType.galaxy_tool if uri.startswith(GALAXY_TOOLS_PREFIX) else RunnableType.galaxy_workflow
    runnable = Runnable(uri, runnable_type)
    return runnable


def cases(runnable: Runnable) -> List["AbstractTestCase"]:
    """Build a `list` of :class:`TestCase` objects for specified runnable."""
    cases: List["AbstractTestCase"] = []

    tests_path = _tests_path(runnable)
    if tests_path is None:
        if runnable.type in (RunnableType.galaxy_tool, RunnableType.galaxy_datamanager):
            if runnable.uri.startswith(GALAXY_TOOLS_PREFIX):
                return [DelayedGalaxyToolTestCase(runnable)]
            tool_source = get_tool_source(runnable.path)
            test_dicts = tool_source.parse_tests_to_dict()
            tool_id = tool_source.parse_id()
            tool_version = tool_source.parse_version()
            for i, test_dict in enumerate(test_dicts.get("tests", [])):
                cases.append(ExternalGalaxyToolTestCase(runnable, tool_id, tool_version, i, test_dict))
        return cases

    return definition_to_test_case(tests_path=tests_path, runnable=runnable)


def definition_to_test_case(tests_path: str, runnable: Runnable) -> List["AbstractTestCase"]:
    with open(tests_path) as f:
        tests_def = yaml.safe_load(f)
    tests_directory = os.path.abspath(os.path.dirname(tests_path))

    def normalize_to_tests_path(path: str) -> str:
        if not os.path.isabs(path):
            absolute_path = os.path.join(tests_directory, path)
        else:
            absolute_path = path
        return os.path.normpath(absolute_path)

    if not isinstance(tests_def, list):
        message = TEST_FILE_NOT_LIST_MESSAGE % tests_path
        raise Exception(message)

    cases: List["AbstractTestCase"] = []
    for i, test_def in enumerate(tests_def):
        if "job" not in test_def:
            message = TEST_FIELD_MISSING_MESSAGE % (i + 1, tests_path, "job")
            raise Exception(message)
        job_def = test_def["job"]
        if isinstance(job_def, dict):
            job_path = None
            job = job_def
        else:
            job_path = normalize_to_tests_path(job_def)
            job = None

        doc = test_def.get("doc")
        output_expectations = test_def.get("outputs", {})
        case = TestCase(
            runnable=runnable,
            tests_directory=tests_directory,
            output_expectations=output_expectations,
            index=i,
            job_path=job_path,
            job=job,
            doc=doc,
        )
        cases.append(case)
    return cases


class AbstractTestCase(metaclass=abc.ABCMeta):
    """Description of a test case for a runnable."""

    def structured_test_data(self, run_response):
        """Result of executing this test case - a "structured_data" dict.

        :rtype: dict
        :return:
                 For example::

                   {
                       "id": "",
                       "has_data": true,
                       "data": {
                           "status": "success", // error, skip,
                           "job": {
                               "command_line": "cat moo",
                               "stdout": "",
                               "stderr": ""
                           },
                           "output_problems": [],
                           "execution_problem": "",
                           "inputs" = {},
                           "problem_log": ""
                       }
                   }
        """


class TestCase(AbstractTestCase):
    """Describe an abstract test case for a specified runnable."""

    def __init__(
        self,
        runnable: Runnable,
        tests_directory: str,
        output_expectations: Dict[str, Any],
        job_path: Optional[str],
        job: Optional[Dict],
        index: int,
        doc: Optional[str],
    ) -> None:
        """Construct TestCase object from required attributes."""
        self.runnable = runnable
        self.job_path = job_path
        self.job = job
        self.output_expectations = output_expectations
        self.tests_directory = tests_directory
        self.index = index
        self.doc = doc

    def __repr__(self) -> str:
        return f"TestCase ({self.doc}) for runnable ({self.runnable}) with job ({self.job}) and expected outputs ({self.output_expectations}) in directory ({self.tests_directory}) with id ({self.index})"

    def structured_test_data(self, run_response: "RunResponse") -> Dict[str, Any]:
        """Check a test case against outputs dictionary."""
        return run_response.structured_data(self)

    @property
    def _job(self):
        if self.job_path is not None:
            with open(self.job_path) as f:
                return yaml.safe_load(f)
        else:
            return self.job

    @property
    def input_ids(self) -> List[str]:
        """Labels of inputs specified in test description."""
        return list(self._job.keys())

    @property
    def tested_output_ids(self) -> List[str]:
        """Labels of outputs checked in test description."""
        return list(self.output_expectations.keys())

    def _check_output(
        self,
        output_id: str,
        output_value: Any,
        output_test: Any,
    ) -> List[str]:
        output_problems = []
        if not isinstance(output_test, dict):
            if output_test != output_value:
                message = f"Output [{output_id}] value [{output_value}] does not match expected value [{output_test}]."
                output_problems.append(message)
        else:
            if not for_collections(output_test):
                if not isinstance(output_value, dict):
                    message = f"Expected file properties for output [{output_id}]"
                    print(message)
                    print(output_value)
                    output_problems.append(message)
                    return output_problems
                if "path" not in output_value and "location" in output_value:
                    assert output_value["location"].startswith("file://")
                    output_value["path"] = output_value["location"][len("file://") :]
                if "path" not in output_value:
                    message = f"No path specified for expected output file [{output_id}]"
                    output_problems.append(message)
                    print(message)
                    return output_problems
            else:
                output_test["name"] = output_id

            output_problems.extend(
                check_output(
                    self.runnable,
                    output_value,
                    output_test,
                    # TODO: needs kwds in here...
                )
            )

        return output_problems

    @property
    def _test_id(self) -> str:
        if self.runnable.type in [
            RunnableType.cwl_tool,
            RunnableType.galaxy_tool,
        ]:
            return get_tool_source(self.runnable.path).parse_id()
        else:
            return os.path.basename(self.runnable.uri)


class ExternalGalaxyToolTestCase(AbstractTestCase):
    """Special class of AbstractCase that doesn't use job_path but uses test data from a Galaxy server."""

    def __init__(
        self,
        runnable: Runnable,
        tool_id: Optional[str],
        tool_version: Optional[str],
        test_index: Optional[int],
        test_dict: Any,
    ) -> None:
        """Construct TestCase object from required attributes."""
        self.runnable = runnable
        self.tool_id = tool_id
        self.tool_version = tool_version
        self.test_index = test_index
        self.test_dict = test_dict

    def structured_test_data(self, run_response: Dict[str, Any]) -> Dict[str, Any]:
        """Just return the structured_test_data generated from galaxy-tool-util for this test variant."""
        return run_response


class DelayedGalaxyToolTestCase(ExternalGalaxyToolTestCase):
    """Special class that requires installing tools prior to finding test cases."""

    def __init__(self, runnable: Runnable) -> None:
        super().__init__(runnable, tool_id=None, tool_version=None, test_index=None, test_dict=None)


def _tests_path(runnable: Runnable) -> Optional[str]:
    if not runnable.is_single_artifact:
        raise NotImplementedError("Tests for directories are not yet implemented.")

    runnable_path = runnable.path
    base, _ = os.path.splitext(runnable_path)

    for test_suffix in TEST_SUFFIXES:
        for test_extension in TEST_EXTENSIONS:
            test_path = base + test_suffix + test_extension
            if os.path.exists(test_path):
                return test_path

    return None


def get_outputs(runnable: Runnable, gi: Optional[GalaxyInstance] = None) -> List["RunnableOutput"]:
    """Return a list of :class:`RunnableOutput` objects for this runnable.

    Supply bioblend user Galaxy instance object (as gi) if additional context
    needed to resolve workflow details.
    """
    if not runnable.is_single_artifact:
        raise NotImplementedError("Cannot generate outputs for a directory.")
    if runnable.type in [RunnableType.galaxy_tool, RunnableType.cwl_tool]:
        tool_source = get_tool_source(runnable.path)
        # TODO: do something with collections at some point
        output_datasets, _ = tool_source.parse_outputs(None)
        return [ToolOutput(o) for o in output_datasets.values()]
    elif runnable.type == RunnableType.galaxy_workflow:
        workflow_outputs = describe_outputs(runnable, gi=gi)
        return [GalaxyWorkflowOutput(o) for o in workflow_outputs]
    elif runnable.type == RunnableType.cwl_workflow:
        workflow = workflow_proxy(runnable.path, strict_cwl_validation=False)
        return [CwlWorkflowOutput(label) for label in workflow.output_labels]
    else:
        raise NotImplementedError("Getting outputs for this artifact type is not yet supported.")


class RunnableOutput(metaclass=abc.ABCMeta):
    """Description of a single output of an execution of a Runnable."""

    @abc.abstractproperty
    def get_id(self):
        """An identifier that describes this output."""

    def is_optional(self):
        return False


class ToolOutput(RunnableOutput):
    """Implementation of RunnableOutput corresponding to Galaxy tool outputs."""

    def __init__(self, tool_output):
        self._tool_output = tool_output

    def get_id(self):
        return self._tool_output.name


class GalaxyWorkflowOutput(RunnableOutput):
    """Implementation of RunnableOutput corresponding to Galaxy workflow outputs."""

    def __init__(self, workflow_output: WorkflowOutput) -> None:
        self._workflow_output = workflow_output

    def get_id(self) -> Optional[str]:
        return self._workflow_output.label

    def is_optional(self):
        return self.workflow_output.optional

    @property
    def workflow_output(self):
        return self._workflow_output


class CwlWorkflowOutput(RunnableOutput):
    """Implementation of RunnableOutput corresponding to CWL outputs."""

    def __init__(self, label: str) -> None:
        self._label = label

    def get_id(self) -> str:
        return self._label


class RunResponse(metaclass=abc.ABCMeta):
    """Description of an attempt for an engine to execute a Runnable."""

    @property
    def start_datetime(self) -> None:
        """Start datetime of run."""
        return None

    @property
    def end_datetime(self) -> None:
        """End datetime of run."""
        return None

    @abc.abstractproperty
    def was_successful(self) -> bool:
        """Indicate whether an error was encountered while executing this runnable.

        If successful, response should conform to the SuccessfulRunResponse interface,
        otherwise it will conform to the ErrorRunResponse interface.
        """

    @abc.abstractproperty
    def job_info(self):
        """If job information is available, return as dictionary."""

    @abc.abstractproperty
    def invocation_details(self):
        """If workflow invocation details are available, return as dictionary."""

    @abc.abstractproperty
    def log(self):
        """If engine related log is available, return as text data."""

    @abc.abstractproperty
    def outputs_dict(self):
        """Return a dict of output descriptions."""

    def get_output(self, output_id):
        """Fetch output from engine."""
        return self.outputs_dict.get(output_id)

    def structured_data(self, test_case: Optional[TestCase] = None) -> Dict[str, Any]:
        output_problems = []
        if self.was_successful:
            execution_problem = None
            if test_case:
                for output_id, output_test in test_case.output_expectations.items():
                    output_value = self.get_output(output_id)
                    if not output_value:
                        message = f"Expected output [{output_id}] not found in results."
                        output_problems.append(message)
                        continue

                    output_problems.extend(test_case._check_output(output_id, output_value, output_test))
            if output_problems:
                status = "failure"
            else:
                status = "success"
        else:
            execution_problem = getattr(self, "error_message", None)
            status = "error"
        data_dict: Dict[str, Any] = dict(status=status)
        if status != "success":
            data_dict["output_problems"] = output_problems
            data_dict["execution_problem"] = execution_problem
        log = self.log
        if log is not None:
            data_dict["problem_log"] = log
        job_info = self.job_info
        if job_info is not None:
            data_dict["job"] = job_info
        invocation_details = self.invocation_details
        if invocation_details is not None:
            data_dict["invocation_details"] = invocation_details
        if self.start_datetime is not None:
            data_dict["start_datetime"] = self.start_datetime.isoformat()
        if self.end_datetime is not None:
            data_dict["end_datetime"] = self.end_datetime.isoformat()
        if test_case:
            data_dict["inputs"] = test_case._job
            return dict(
                id=(f"{test_case._test_id}_{test_case.index}"),
                has_data=True,
                data=data_dict,
                doc=test_case.doc,
                test_type=test_case.runnable.type.name,
            )
        else:
            assert isinstance(self, SuccessfulRunResponse)
            return dict(
                id=self._runnable.uri,
                has_data=True,
                data=data_dict,
                doc=None,
                test_type=self._runnable.type.name,
            )


class SuccessfulRunResponse(RunResponse, metaclass=abc.ABCMeta):
    """Description of the results of an engine executing a Runnable."""

    def __init__(self, runnable: "Runnable") -> None:
        self._runnable = runnable

    @property
    def was_successful(self):
        """Return `True` to indicate this run was successful."""
        return True


class ErrorRunResponse(RunResponse):
    """Description of an error while attempting to execute a Runnable."""

    def __init__(
        self, error_message, job_info=None, invocation_details=None, log=None, start_datetime=None, end_datetime=None
    ):
        """Create an ErrorRunResponse with specified error message."""
        self._error_message = error_message
        self._job_info = job_info
        self._invocation_details = invocation_details
        self._log = log
        self._start_datetime = start_datetime
        self._end_datetime = end_datetime

    @property
    def start_datetime(self):
        """Start datetime of run."""
        return self._start_datetime

    @property
    def end_datetime(self):
        """End datetime of run."""
        return self._end_datetime

    @property
    def error_message(self):
        """Error message describing the problem with execution of the runnable."""
        return self._error_message

    @property
    def was_successful(self):
        """Return `False` to indicate this run was successful."""
        return False

    @property
    def job_info(self):
        """Return potentially null stored `job_info` dict."""
        return self._job_info

    @property
    def invocation_details(self):
        return self._invocation_details

    @property
    def log(self):
        """Return potentially null stored `log` text."""
        return self._log

    @property
    def outputs_dict(self):
        return {}

    def __str__(self):
        """Print a helpful error description of run."""
        message = f"Run failed with message [{self.error_message}]"
        log = self.log
        if log:
            message += f" and log [{log}]"
        return message


__all__ = (
    "cases",
    "ErrorRunResponse",
    "for_path",
    "for_paths",
    "get_outputs",
    "Runnable",
    "RunnableType",
    "RunResponse",
    "RunnableOutput",
    "SuccessfulRunResponse",
    "TestCase",
)
