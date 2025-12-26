"""Workflow linting entry point - main script."""
import argparse
import os
import sys

from gxformat2._scripts import ensure_format2
from gxformat2.linting import LintContext
from gxformat2.markdown_parse import validate_galaxy_markdown
from gxformat2.normalize import Inputs
from gxformat2.yaml import ordered_load, ordered_load_path

EXIT_CODE_SUCCESS = 0
EXIT_CODE_LINT_FAILED = 1
EXIT_CODE_FORMAT_ERROR = 2
EXIT_CODE_FILE_PARSE_FAILED = 3

LINT_FAILED_NO_OUTPUTS = "Workflow contained no outputs"
LINT_FAILED_OUTPUT_NO_LABEL = "Workflow contained output without a label"


def ensure_key(lint_context, has_keys, key, has_class=None, has_value=None):
    if key not in has_keys:
        lint_context.error("expected to find key [{key}] but absent", key=key)
        return None

    value = has_keys[key]
    return ensure_key_has_value(lint_context, has_keys, key, value, has_class=has_class, has_value=has_value)


def ensure_key_if_present(lint_context, has_keys, key, default=None, has_class=None):
    if key not in has_keys:
        return default

    value = has_keys[key]
    return ensure_key_has_value(lint_context, has_keys, key, value, has_class=has_class, has_value=None)


def ensure_key_has_value(lint_context, has_keys, key, value, has_class=None, has_value=None):
    if has_class is not None and not isinstance(value, has_class):
        lint_context.error(f"expected value [{value}] with key [{key}] to be of class {has_class}")
    if has_value is not None and value != has_value:
        lint_context.error(f"expected value [{value}] with key [{key}] to be {has_value}")
    return value


def _lint_step_errors(lint_context, step):
    step_errors = step.get("errors")
    if step_errors is not None:
        lint_context.warn(f"tool step contains error indicated during Galaxy export - {step_errors}")


def lint_ga_path(lint_context, path):
    """Apply linting of native workflows to specified path."""
    workflow_dict = ordered_load_path(path)
    return lint_ga(lint_context, workflow_dict, path=path)


def lint_ga(lint_context, workflow_dict, path=None):
    """Lint a native/legacy style Galaxy workflow and populate the corresponding LintContext."""
    ensure_key(lint_context, workflow_dict, "format-version", has_value="0.1")
    ensure_key(lint_context, workflow_dict, "a_galaxy_workflow", has_value="true")

    native_steps = ensure_key(lint_context, workflow_dict, "steps", has_class=dict) or {}

    found_outputs = False
    found_output_without_label = False
    for order_index_str, step in native_steps.items():
        if not order_index_str.isdigit():
            lint_context.error("expected step_key to be integer not [{value}]", value=order_index_str)

        workflow_outputs = ensure_key_if_present(lint_context, step, "workflow_outputs", default=[], has_class=list)
        for workflow_output in workflow_outputs:
            found_outputs = True

            if not workflow_output.get("label"):
                found_output_without_label = True

        step_type = step.get("type")
        if step_type == "subworkflow":
            subworkflow = ensure_key(lint_context, step, "subworkflow", has_class=dict)
            lint_ga(lint_context, subworkflow)

        _lint_step_errors(lint_context, step)
        _lint_tool_if_present(lint_context, step)

    _validate_report(lint_context, workflow_dict)
    if not found_outputs:
        lint_context.warn(LINT_FAILED_NO_OUTPUTS)

    if found_output_without_label:
        lint_context.warn(LINT_FAILED_OUTPUT_NO_LABEL)

    _lint_training(lint_context, workflow_dict)


def lint_format2(lint_context, workflow_dict, path=None):
    """Lint a Format 2 Galaxy workflow and populate the corresponding LintContext."""
    from gxformat2.schema.v19_09 import load_document
    from schema_salad.exceptions import SchemaSaladException  # type: ignore
    try:
        load_document("file://" + os.path.abspath(path))
    except SchemaSaladException as e:
        lint_context.error("Validation failed " + str(e))

    steps = ensure_key_if_present(lint_context, workflow_dict, 'steps', default={}, has_class=(dict, list))
    steps = steps.values() if isinstance(steps, dict) else steps
    for step in steps:
        _lint_step_errors(lint_context, step)
        _lint_tool_if_present(lint_context, step)

    _validate_input_types(lint_context, workflow_dict)
    _validate_report(lint_context, workflow_dict)
    _lint_training(lint_context, workflow_dict)


def _validate_input_types(lint_context: LintContext, workflow_dict: dict):
    try:
        inputs = Inputs(workflow_dict)
    except Exception:
        # bad document, can't process inputs...
        return
    for input_def in inputs._inputs:
        input_type = input_def.get("type")
        if "default" in input_def:
            input_default = input_def['default']
            if input_type == "int":
                if not isinstance(input_default, int):
                    lint_context.error('Input default is of invalid type')
            elif input_type == "float":
                if not isinstance(input_default, (int, float)):
                    lint_context.error('Input default is of invalid type')
            elif input_type == "string":
                if not isinstance(input_default, str):
                    lint_context.error('Input default is of invalid type')


def _lint_tool_if_present(lint_context, step_dict):
    tool_id = step_dict.get('tool_id')
    if tool_id and 'testtoolshed' in tool_id:
        lint_context.warn('Step references a tool from the test tool shed, this should be replaced with a production tool')


def _validate_report(lint_context, workflow_dict):
    report_dict = ensure_key_if_present(lint_context, workflow_dict, "report", default=None, has_class=dict)
    if report_dict is not None:
        markdown = ensure_key(lint_context, report_dict, "markdown", has_class=str)
        if isinstance(markdown, str):
            try:
                validate_galaxy_markdown(markdown)
            except ValueError as e:
                lint_context.error(f"Report markdown validation failed [{e}]")


def _lint_training(lint_context, workflow_dict):
    if lint_context.training_topic is None:
        return

    if "tags" not in workflow_dict:
        lint_context.warn("Missing tag(s).")
    else:
        tags = workflow_dict["tags"]
        if lint_context.training_topic not in tags:
            lint_context.warn(f"Missing expected training topic ({lint_context.training_topic}) as workflow tag.")
    # Move up into individual lints - all workflows should have docs.
    format2_dict = ensure_format2(workflow_dict)
    if "doc" not in format2_dict:
        lint_context.warn("Missing workflow documentation (annotation or doc element)")
    elif not format2_dict["doc"]:
        lint_context.warn("Empty workflow documentation (annotation or doc element)")


def main(argv=None):
    """Script entry point for linting workflows."""
    if argv is None:
        argv = sys.argv
    args = _parser().parse_args(argv[1:])
    path = args.path
    with open(path) as f:
        try:
            workflow_dict = ordered_load(f)
        except Exception:
            return EXIT_CODE_FILE_PARSE_FAILED
    workflow_class = workflow_dict.get("class")
    lint_func = lint_format2 if workflow_class == "GalaxyWorkflow" else lint_ga
    lint_context = LintContext(training_topic=args.training_topic)
    lint_func(lint_context, workflow_dict, path=path)
    lint_context.print_messages()
    if lint_context.found_errors:
        return EXIT_CODE_FORMAT_ERROR
    elif lint_context.found_warns:
        return EXIT_CODE_LINT_FAILED
    else:
        return EXIT_CODE_SUCCESS


def _parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("--training-topic",
                        required=False,
                        help='If this is a training workflow, specify a training topic.')
    parser.add_argument('path', metavar='PATH', type=str,
                        help='workflow path')
    return parser


if __name__ == "__main__":
    sys.exit(main())


__all__ = ('main', 'lint_format2', 'lint_ga')
