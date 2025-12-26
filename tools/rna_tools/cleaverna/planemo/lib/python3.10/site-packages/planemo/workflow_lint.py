import inspect
import json
import os
import re
from collections import OrderedDict
from typing import (
    Any,
    Dict,
    Iterable,
    Iterator,
    List,
    Optional,
    Tuple,
    TYPE_CHECKING,
)
from urllib.parse import urlparse

import requests
import yaml
from bioblend import toolshed
from bioblend.toolshed import ToolShedInstance
from galaxy.tool_util.lint import LintContext
from galaxy.tool_util.loader_directory import EXCLUDE_WALK_DIRS
from galaxy.tool_util.parser.yaml import __to_test_assert_list
from galaxy.tool_util.verify import asserts
from gxformat2.lint import (
    lint_format2,
    lint_ga,
)
from gxformat2.yaml import ordered_load

from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK,
)
from planemo.galaxy.workflows import (
    input_labels,
    MAIN_TOOLSHED_URL,
    output_labels,
    required_input_labels,
)
from planemo.lint import build_lint_args
from planemo.runnable import (
    cases,
    for_path,
    TestCase,
)
from planemo.shed import DOCKSTORE_REGISTRY_CONF

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

POTENTIAL_WORKFLOW_FILES = re.compile(r"^.*(\.yml|\.yaml|\.ga)$")
DOCKSTORE_REGISTRY_CONF_VERSION = "1.2"


class WorkflowLintContext(LintContext):
    # Setup training topic for linting - probably should pass this through
    # from click arguments.
    training_topic = None


def build_wf_lint_args(ctx: "PlanemoCliContext", **kwds) -> Dict[str, Any]:
    lint_args = build_lint_args(ctx, **kwds)
    lint_args["iwc_grade"] = kwds.get("iwc", False)
    return lint_args


def generate_dockstore_yaml(directory: str, publish: bool = True) -> str:
    workflows = []
    all_workflow_paths = list(find_workflow_descriptions(directory))
    for workflow_path in all_workflow_paths:
        test_parameter_path = f"{workflow_path.rsplit('.', 1)[0]}-tests.yml"
        workflow_entry: Dict[str, Any] = {
            # TODO: support CWL
            "name": "main" if len(all_workflow_paths) == 1 else os.path.basename(workflow_path).split(".ga")[0],
            "subclass": "Galaxy",
            "publish": publish,
            "primaryDescriptorPath": f"/{os.path.relpath(workflow_path, directory)}",
        }
        if os.path.exists(test_parameter_path):
            workflow_entry["testParameterFiles"] = [f"/{os.path.relpath(test_parameter_path, directory)}"]

        with open(workflow_path) as f:
            workflow_dict = ordered_load(f)

        # add author
        if len(workflow_dict.get("creator", [])) > 0:
            creators = workflow_dict["creator"]
            authors = []
            for creator in creators:
                author = {}
                # Put name first
                if "name" in creator:
                    author["name"] = creator["name"]
                for field, value in creator.items():
                    if field in ["class", "name"]:
                        continue
                    if field == "identifier":
                        # Check if it is an orcid:
                        # Read https://support.orcid.org/hc/en-us/articles/360006897674-Structure-of-the-ORCID-Identifier
                        orcid = re.findall(r"(?:\d{4}-){3}\d{3}[0-9X]", value)
                        if len(orcid) > 0:
                            # Check the orcid is valid
                            if (
                                requests.get(
                                    f"https://orcid.org/{orcid[0]}",
                                    headers={"Accept": "application/xml"},
                                    allow_redirects=True,
                                ).status_code
                                == 200
                            ):
                                author["orcid"] = orcid[0]
                            else:
                                # If it is not put as identifier
                                author["identifier"] = value
                        else:
                            author["identifier"] = value
                    else:
                        author[field] = value
                authors.append(author)
            workflow_entry["authors"] = authors

        workflows.append(workflow_entry)
    # Force version to the top of file but serializing rest of config separately
    contents = "version: %s\n" % DOCKSTORE_REGISTRY_CONF_VERSION
    contents += yaml.dump({"workflows": workflows}, sort_keys=False)
    return contents


def lint_workflow_artifacts_on_paths(ctx: "PlanemoCliContext", paths: Iterable[str], lint_args: Dict[str, Any]) -> int:
    report_level = lint_args["level"]
    lint_context = WorkflowLintContext(report_level, skip_types=lint_args["skip_types"])
    for path in paths:
        _lint_workflow_artifacts_on_path(lint_context, path, lint_args)

    if lint_context.failed(lint_args["fail_level"]):
        return EXIT_CODE_GENERIC_FAILURE
    else:
        return EXIT_CODE_OK


def _lint_workflow_artifacts_on_path(lint_context: WorkflowLintContext, path: str, lint_args: Dict[str, Any]) -> None:
    if lint_args["iwc_grade"]:
        if not os.path.isdir(path):
            path = os.path.dirname(path)
        lint_context.lint("lint_required_files", _lint_required_files_workflow_dir, path)
        lint_context.lint("lint_changelog", _lint_changelog_version, path)

    for potential_workflow_artifact_path in find_potential_workflow_files(path):
        if os.path.basename(potential_workflow_artifact_path) == DOCKSTORE_REGISTRY_CONF:
            lint_context.lint("lint_dockstore", _lint_dockstore_config, potential_workflow_artifact_path)
            if lint_args["iwc_grade"]:
                lint_context.lint(
                    "lint_dockstore_best_practices",
                    _lint_dockstore_config_best_practices,
                    potential_workflow_artifact_path,
                )

        elif looks_like_a_workflow(potential_workflow_artifact_path):

            def structure(path, lint_context):
                with open(path) as f:
                    workflow_dict = ordered_load(f)
                workflow_class = workflow_dict.get("class")
                lint_func = lint_format2 if workflow_class == "GalaxyWorkflow" else lint_ga
                lint_func(lint_context, workflow_dict, path=path)

            lint_context.lint("lint_structure", structure, potential_workflow_artifact_path)
            if lint_args["iwc_grade"]:
                lint_context.lint("lint_release", _lint_release, potential_workflow_artifact_path)
            lint_context.lint("lint_best_practices", _lint_best_practices, potential_workflow_artifact_path)
            lint_context.lint("lint_tests", _lint_tsts, potential_workflow_artifact_path)
            lint_context.lint("lint_tool_ids", _lint_tool_ids, potential_workflow_artifact_path)
        else:
            # Allow linting ro crates and such also
            pass


# misspell for pytest
def _lint_tsts(path: str, lint_context: WorkflowLintContext) -> None:
    runnables = for_path(path)
    if not isinstance(runnables, list):
        runnables = [runnables]
    for runnable in runnables:
        test_cases = cases(runnable)
        if len(test_cases) == 0:
            lint_context.warn("Workflow missing test cases.")
            return
        all_tests_valid = True
        for test_case in test_cases:
            if isinstance(test_case, TestCase):
                if not _lint_case(path, test_case, lint_context):
                    all_tests_valid = False
            else:
                lint_context.warn(f"Test case of type {type(test_case)} not currently supported.")
                all_tests_valid = False

        if all_tests_valid:
            lint_context.valid(f"Tests appear structurally correct for {runnable.path}")


def _lint_best_practices(path: str, lint_context: WorkflowLintContext) -> None:  # noqa: C901
    """
    This function duplicates the checks made by Galaxy's best practices panel:
    https://github.com/galaxyproject/galaxy/blob/5396bb15fe8cfcf2e89d46c1d061c49b60e2f0b1/client/src/components/Workflow/Editor/Lint.vue
    """

    def check_json_for_untyped_params(j):
        values = j.values() if isinstance(j, dict) else j
        for value in values:
            if type(value) in [list, dict, OrderedDict]:
                if check_json_for_untyped_params(value):
                    return True
            elif isinstance(value, str):
                if re.match(r"\$\{.+?\}", value):
                    return True
        return False

    with open(path) as f:
        workflow_dict = ordered_load(f)

    steps = workflow_dict.get("steps", {})

    # annotation
    if not workflow_dict.get("annotation"):
        lint_context.warn("Workflow is not annotated.")

    # creator
    creators = workflow_dict.get("creator", [])
    if not len(creators) > 0:
        lint_context.warn("Workflow does not specify a creator.")
    else:
        if not isinstance(creators, list):
            # Don't know if this can happen, if we implement schema validation on the Galaxy side
            # this won't be needed.
            creators = [creators]
        for creator in creators:
            if creator.get("class", "").lower() == "person" and "identifier" in creator:
                identifier = creator["identifier"]
                parsed_url = urlparse(identifier)
                if not parsed_url.scheme:
                    lint_context.warn(
                        f'Creator identifier "{identifier}" should be a fully qualified URI, for example "https://orcid.org/0000-0002-1825-0097".'
                    )

    # license
    if not workflow_dict.get("license"):
        lint_context.warn("Workflow does not specify a license.")

    # checks on individual steps
    for step in steps.values():
        # disconnected inputs
        if step.get("type") not in ["data_collection_input", "parameter_input"]:
            for input in step.get("inputs", []):
                if input.get("name") not in step.get("input_connections"):  # TODO: check optional
                    lint_context.warn(
                        f"Input {input.get('name')} of workflow step {step.get('annotation') or step.get('id')} is disconnected."
                    )

        # missing metadata
        if not step.get("annotation"):
            lint_context.warn(f"Workflow step with ID {step.get('id')} has no annotation.")
        if not step.get("label"):
            lint_context.warn(f"Workflow step with ID {step.get('id')} has no label.")

        # untyped parameters
        if workflow_dict.get("class") == "GalaxyWorkflow":
            tool_state = step.get("tool_state", {})
            pjas = step.get("out", {})
        else:
            tool_state = json.loads(step.get("tool_state", "{}"))
            pjas = step.get("post_job_actions", {})

        if check_json_for_untyped_params(tool_state):
            lint_context.warn(f"Workflow step with ID {step.get('id')} specifies an untyped parameter as an input.")

        if check_json_for_untyped_params(pjas):
            lint_context.warn(
                f"Workflow step with ID {step.get('id')} specifies an untyped parameter in the post-job actions."
            )

        # unlabeled outputs are checked by gxformat2, no need to check here


def _lint_case(path: str, test_case: TestCase, lint_context: WorkflowLintContext) -> bool:
    test_valid = True

    i_labels = input_labels(workflow_path=path)
    job_keys = test_case.input_ids
    for key in job_keys:
        if key not in i_labels:
            # consider an error instead?
            lint_context.warn(
                f"Unknown workflow input in test job definition [{key}], workflow inputs are [{i_labels}]"
            )
            test_valid = False

    # check non-optional parameters are set
    for required_label in required_input_labels(path):
        if required_label not in job_keys:
            template = "Non-optional input has no value specified in workflow test job [%s], job specifies inputs [%s]"
            lint_context.error(template % (required_label, job_keys))
            test_valid = False

    for input_id, input_def in test_case._job.items():
        if not _tst_input_valid(test_case, input_id, input_def, lint_context):
            test_valid = False

    test_output_ids = test_case.tested_output_ids
    o_labels = output_labels(path)
    found_valid_expectation = False
    for test_output_id in test_output_ids:
        if test_output_id not in o_labels:
            template = "Test found for unknown workflow output [%s], workflow outputs [%s]"
            lint_context.error(template % (test_output_id, o_labels))
            test_valid = False
        else:
            found_valid_expectation = True
        # TODO: validate structure of test expectations

        output_expectations = test_case.output_expectations[test_output_id]
        test_valid = is_valid_output_expectations(lint_context=lint_context, output_expectations=output_expectations)

    if not found_valid_expectation:
        lint_context.warn("Found no valid test expectations for workflow test")
        test_valid = False

    return test_valid


def is_valid_output_expectations(lint_context, output_expectations):
    all_assertion_definitions = []
    if isinstance(output_expectations, (int, str, float, bool)):
        # CWL style parameter output
        return True
    elif "element_tests" in output_expectations:
        # This is a collection
        for element_id in output_expectations["element_tests"]:
            all_assertion_definitions.append(output_expectations["element_tests"][element_id].get("asserts"))
    else:
        all_assertion_definitions.append(output_expectations.get("asserts"))
    for assertion_definitions in all_assertion_definitions:
        if not _check_test_assertions(lint_context, assertion_definitions):
            return False
    return True


def _check_test_assertions(
    lint_context: WorkflowLintContext, assertion_definitions: Optional[Dict[str, Dict[str, Any]]]
) -> bool:
    # we are already in Python, not XML, so it is simpler to lint assertions by checking against the
    # Python functions directly, rather than checking against galaxy.xsd as for tool linting
    assertions_valid = True
    if assertion_definitions:
        # Can be either a dictionary with different assert
        # Or a list with max of asserts and potentially identical
        # We transform to list:
        assertion_definitions_list = __to_test_assert_list(assertion_definitions)
        for assertion_description in assertion_definitions_list:
            function = asserts.assertion_functions.get(f"assert_{assertion_description['tag']}")
            if function is None:
                lint_context.error(f"Invalid assertion: assert_{assertion_description['tag']} does not exists")
                assertions_valid = False
                continue
            signature = inspect.signature(function)
            function_args = inspect.getfullargspec(function).args
            assertion_params = assertion_description["attributes"].copy()
            if "verify_assertions_function" in function_args:
                assertion_params["verify_assertions_function"] = asserts.verify_assertion
            if "children" in function_args:
                if "asserts" in assertion_params:
                    assertion_params["children"] = assertion_params["asserts"]
                    del assertion_params["asserts"]
                if "children" not in assertion_params:
                    assertion_params["children"] = []
                _check_test_assertions(lint_context, assertion_params["children"])
            del assertion_params["that"]
            try:
                # try mapping the function with the attributes supplied and check for TypeError
                signature.bind("", **assertion_params)
            except AssertionError:
                pass
            except TypeError as e:
                lint_context.error(f"Invalid assertion in tests: assert_{assertion_description['tag']} {str(e)}")
                assertions_valid = False
    return assertions_valid


def _tst_input_valid(
    test_case: TestCase,
    input_id: str,
    input_def: Dict[str, Any],
    lint_context: WorkflowLintContext,
) -> bool:
    if isinstance(input_def, dict):  # else assume it is a parameter
        clazz = input_def.get("class")
        if clazz == "File":
            input_path = input_def.get("path")
            if input_path:
                if not os.path.isabs(input_path):
                    input_path = os.path.join(test_case.tests_directory, input_path)
                if not os.path.exists(input_path):
                    message = "Test referenced File path [%s] not found" % input_path
                    lint_context.warn(message)
                    return False
        elif clazz == "Collection":
            for elem in input_def.get("elements", []):
                elem_valid = _tst_input_valid(test_case, input_id, elem, lint_context)
                if not elem_valid:
                    return False
    return True


def _lint_dockstore_config(path: str, lint_context: WorkflowLintContext) -> None:
    dockstore_yaml = None
    try:
        with open(path) as f:
            dockstore_yaml = yaml.safe_load(f)
    except Exception:
        lint_context.error("Invalid YAML found in %s" % DOCKSTORE_REGISTRY_CONF)
        return

    if not isinstance(dockstore_yaml, dict):
        lint_context.error("Invalid YAML contents found in %s" % DOCKSTORE_REGISTRY_CONF)
        return

    if "version" not in dockstore_yaml:
        lint_context.error("Invalid YAML contents found in %s, no version defined" % DOCKSTORE_REGISTRY_CONF)
    if str(dockstore_yaml.get("version")) != DOCKSTORE_REGISTRY_CONF_VERSION:
        lint_context.error("Invalid YAML version found in %s." % DOCKSTORE_REGISTRY_CONF)

    if "workflows" not in dockstore_yaml:
        lint_context.error("Invalid YAML contents found in %s, no workflows defined" % DOCKSTORE_REGISTRY_CONF)
        return

    workflow_entries = dockstore_yaml.get("workflows")
    if not isinstance(workflow_entries, list):
        lint_context.error("Invalid YAML contents found in %s, workflows not a list" % DOCKSTORE_REGISTRY_CONF)
        return

    if len(workflow_entries) == 0:
        lint_context.error("No workflow specified in the .dockstore.yml.")

    workflow_names_in_dockstore = []
    for workflow_entry in workflow_entries:
        _lint_dockstore_workflow_entry(lint_context, os.path.dirname(path), workflow_entry)
        workflow_name = workflow_entry.get("name", "")
        if workflow_name in workflow_names_in_dockstore:
            lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} has multiple workflow entries with the same name")
        workflow_names_in_dockstore.append(workflow_name)


def _lint_dockstore_workflow_entry(
    lint_context: WorkflowLintContext, directory: str, workflow_entry: Dict[str, Any]
) -> None:
    if not isinstance(workflow_entry, dict):
        lint_context.error("Invalid YAML contents found in %s, workflow entry not a dict" % DOCKSTORE_REGISTRY_CONF)
        return

    found_errors = False
    for required_key in ["primaryDescriptorPath", "subclass"]:
        if required_key not in workflow_entry:
            lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} workflow entry missing required key {required_key}")
            found_errors = True

    if found_errors:
        # Don't do the rest of the validation for a broken file.
        return

    if workflow_entry.get("subclass") != "Galaxy":
        lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} workflow entry subclass must be 'Galaxy'.")

    descriptor_path = workflow_entry["primaryDescriptorPath"]
    test_files = workflow_entry.get("testParameterFiles", [])

    for referenced_file in [descriptor_path] + test_files:
        referenced_path = os.path.join(directory, referenced_file[1:])
        if not os.path.exists(referenced_path):
            lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} workflow entry references absent file {referenced_file}")

    # Check there is no space in name:
    workflow_name = workflow_entry.get("name", "")
    # Check the name has no space
    if " " in workflow_name:
        lint_context.error(
            "Dockstore does not accept workflow names with space.",
            f"Change '{workflow_name}' in {DOCKSTORE_REGISTRY_CONF}.",
        )

    # Check there is not mailto
    for author in workflow_entry.get("authors", []):
        if author.get("email", "").startswith("mailto:"):
            lint_context.error("email field of the .dockstore.yml must not contain 'mailto:'")


def looks_like_a_workflow(path: str) -> bool:
    """Return boolean indicating if this path looks like a workflow."""
    if POTENTIAL_WORKFLOW_FILES.match(os.path.basename(path)):
        with open(path) as f:
            workflow_dict = ordered_load(f)
        if not isinstance(workflow_dict, dict):
            # Not exactly right - could have a #main def - do better and sync with Galaxy.
            return False
        return bool(workflow_dict.get("class") == "GalaxyWorkflow" or workflow_dict.get("a_galaxy_workflow"))
    return False


def find_workflow_descriptions(directory: str) -> Iterator[str]:
    for potential_workflow_artifact_path in find_potential_workflow_files(directory):
        if looks_like_a_workflow(potential_workflow_artifact_path):
            yield potential_workflow_artifact_path


def find_potential_workflow_files(directory: str) -> List[str]:
    """Return a list of potential workflow files in a directory."""
    if not os.path.exists(directory):
        raise ValueError(f"Directory not found {directory}")

    matches = []
    if os.path.isdir(directory):
        for root, dirnames, filenames in os.walk(directory):
            # exclude some directories (like .hg) from traversing
            dirnames[:] = [dir for dir in dirnames if dir not in EXCLUDE_WALK_DIRS]
            for filename in filenames:
                if POTENTIAL_WORKFLOW_FILES.match(filename):
                    matches.append(os.path.join(root, filename))
    else:
        matches.append(directory)
    return matches


def find_repos_from_tool_id(tool_id: str, ts: ToolShedInstance) -> Tuple[str, Dict[str, Any]]:
    """
    Return a string which indicates what failed and dict with all revisions for a given tool id
    """
    if "/repos" not in tool_id:
        return ("", {})  # assume a built in tool
    *_, owner, name, _tool_id, _version = tool_id.split("/")

    try:
        repo = ts.repositories.get_repositories(name, owner)[0]
        repos = ts.repositories._get(url=f"{ts.repositories._make_url()}/{repo['id']}/metadata")
    except Exception as e:
        return (f"The ToolShed returned an error when searching for the most recent version of {tool_id}: {e}", {})
    if len(repos) == 0:
        return (f"The tool {tool_id} is not in the toolshed (may have been tagged as invalid).", {})
    else:
        return ("", repos)


def assert_valid_tool_id_in_tool_shed(tool_id: str, ts: ToolShedInstance) -> Optional[str]:
    if "/repos" not in tool_id:
        return None
    warning_msg, repos = find_repos_from_tool_id(tool_id, ts)
    if warning_msg:
        return warning_msg
    for repo in repos.values():
        tools = repo.get("tools", [])
        for tool in tools:
            if tool_id == tool.get("guid"):
                return None
    return f"The tool {tool_id} is not in the toolshed (may have been tagged as invalid)."


def _lint_tool_ids(path: str, lint_context: WorkflowLintContext) -> None:
    def _lint_tool_ids_steps(lint_context: WorkflowLintContext, wf_dict: Dict, ts: ToolShedInstance) -> bool:
        """Returns whether a single tool_id was invalid"""
        failed = False
        steps = wf_dict.get("steps", {})
        for step in steps.values():
            if step.get("type", "tool") == "tool" and not step.get("run", {}).get("class") == "GalaxyWorkflow":
                warning_msg = assert_valid_tool_id_in_tool_shed(step["tool_id"], ts)
                if warning_msg:
                    lint_context.error(warning_msg)
                    failed = True
            elif step.get("type") == "subworkflow":  # GA SWF
                sub_failed = _lint_tool_ids_steps(lint_context, step["subworkflow"], ts)
                if sub_failed:
                    failed = True
            elif step.get("run", {}).get("class") == "GalaxyWorkflow":  # gxformat2 SWF
                sub_failed = _lint_tool_ids_steps(lint_context, step["run"], ts)
                if sub_failed:
                    failed = True
            else:
                continue
        return failed

    with open(path) as f:
        workflow_dict = ordered_load(f)
    ts = toolshed.ToolShedInstance(url=MAIN_TOOLSHED_URL)
    failed = _lint_tool_ids_steps(lint_context, workflow_dict, ts)
    if not failed:
        lint_context.valid("All tool ids appear to be valid.")
    return None


def _lint_required_files_workflow_dir(path: str, lint_context: WorkflowLintContext) -> None:
    # Check all required files are present
    required_files = ["README.md", "CHANGELOG.md", ".dockstore.yml"]
    for required_file in required_files:
        if not os.path.exists(os.path.join(path, required_file)):
            lint_context.error(f"The file {required_file} is missing but required.")


def _get_changelog_version(path: str) -> str:
    # Get the version from the CHANGELOG.md
    version = ""
    if not os.path.exists(os.path.join(path, "CHANGELOG.md")):
        return version
    with open(os.path.join(path, "CHANGELOG.md"), "r") as f:
        for line in f:
            if line.startswith("## ["):
                version = line.split("]")[0].replace("## [", "")
                break
    return version


def _lint_changelog_version(path: str, lint_context: WorkflowLintContext) -> None:
    # Check the version can be get from the CHANGELOG.md
    if not os.path.exists(os.path.join(path, "CHANGELOG.md")):
        return
    if _get_changelog_version(path) == "":
        lint_context.error(
            "No version found in CHANGELOG. The version should be in a line that starts like '## [version number]'"
        )


def _lint_release(path, lint_context):
    with open(path) as f:
        workflow_dict = ordered_load(f)
    if "release" not in workflow_dict:
        lint_context.error(f"The workflow {path} has no release")
    else:
        # Try to get the version from the CHANGELOG:
        version = _get_changelog_version(os.path.dirname(path))
        if version != "" and workflow_dict.get("release") != version:
            lint_context.error(f"The release of workflow {path} does not match the version in the CHANGELOG.")


def _lint_dockstore_config_best_practices(path: str, lint_context: WorkflowLintContext) -> None:
    dockstore_yaml = None
    try:
        with open(path) as f:
            dockstore_yaml = yaml.safe_load(f)
    except Exception:
        return

    if not isinstance(dockstore_yaml, dict):
        return

    workflow_entries = dockstore_yaml.get("workflows")
    if not isinstance(workflow_entries, list):
        return

    for workflow_entry in workflow_entries:
        _lint_dockstore_workflow_entry_best_practices(lint_context, os.path.dirname(path), workflow_entry)


def _lint_dockstore_workflow_entry_best_practices(
    lint_context: WorkflowLintContext, directory: str, workflow_entry: Dict[str, Any]
) -> None:
    for recommended_key in ["testParameterFiles", "name"]:
        if recommended_key not in workflow_entry:
            lint_context.error(f"{DOCKSTORE_REGISTRY_CONF} workflow entry missing recommended key {recommended_key}")

    workflow_name = workflow_entry.get("name", "")
    # Check there is at least one author
    if len(workflow_entry.get("authors", [])) == 0:
        lint_context.error(f"Workflow {workflow_name} have no 'authors' in the {DOCKSTORE_REGISTRY_CONF}.")
