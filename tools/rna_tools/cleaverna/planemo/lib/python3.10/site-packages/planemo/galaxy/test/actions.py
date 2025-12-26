"""Actions related to running and reporting on Galaxy-specific testing."""

import json

import click
from galaxy.util import unicodify

from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_NO_SUCH_TARGET,
    EXIT_CODE_OK,
)
from planemo.io import (
    info,
    warn,
)
from planemo.reports import (
    allure,
    build_report,
)
from planemo.test.results import get_dict_value
from . import structures as test_structures

NO_XUNIT_REPORT_MESSAGE = (
    "Cannot locate xUnit report [%s] for tests - required to build planemo report and summarize tests."
)
NO_JSON_REPORT_MESSAGE = (
    "Cannot locate JSON report [%s] for tests - required to build planemo report and summarize tests."
)
REPORT_NOT_CHANGED = (
    "Galaxy failed to update test report [%s] for tests - required to build planemo report and summarize tests."
)
NO_TESTS_MESSAGE = "No tests were executed - see Galaxy output for details."
ALL_TESTS_PASSED_MESSAGE = "All %d test(s) executed passed."
PROBLEM_COUNT_MESSAGE = (
    "There were problems with %d test(s) - out of %d test(s) executed. See %s for detailed breakdown."
)
GENERIC_PROBLEMS_MESSAGE = "One or more tests failed. See %s for detailed breakdown."
GENERIC_TESTS_PASSED_MESSAGE = "No failing tests encountered."
TEST_DATA_UPDATED_MESSAGE = "Test data were updated and tests were rerun."
TEST_DATA_NOT_UPDATED_MESSAGE = "%s Therefore, no test data were updated." % ALL_TESTS_PASSED_MESSAGE


def handle_reports_and_summary(ctx, structured_data, exit_code=None, kwds=None):
    """Produce reports and print summary, return 0 if tests passed.

    If ``exit_code`` is set - use underlying test source for return
    code and test success determination, otherwise infer from supplied
    test data.
    """
    if kwds is None:
        kwds = {}
    handle_reports(ctx, structured_data, kwds)
    summary_exit_code = _handle_summary(structured_data, **kwds)
    return exit_code if exit_code is not None else summary_exit_code


def merge_reports(input_paths, output_path):
    reports = []
    for path in input_paths:
        with open(path, encoding="utf-8") as f:
            reports.append(json.load(f))
    tests = []
    for report in reports:
        tests.extend(report["tests"])
    tests = sorted(tests, key=lambda k: k["id"])
    merged_report = {"tests": tests}
    with open(output_path, mode="w", encoding="utf-8") as out:
        out.write(unicodify(json.dumps(merged_report)))


def handle_reports(ctx, structured_data, kwds):
    """Write reports based on user specified kwds."""
    exceptions = []
    structured_report_file = kwds.get("test_output_json", None)
    if structured_report_file:
        try:
            with open(structured_report_file, mode="w", encoding="utf-8") as f:
                f.write(unicodify(json.dumps(structured_data, indent=4, sort_keys=True)))
        except Exception as e:
            exceptions.append(e)

    for report_type in ["html", "markdown", "markdown_minimal", "text", "xunit", "junit", "allure"]:
        try:
            _handle_test_output_file(ctx, report_type, structured_data, kwds)
        except Exception as e:
            exceptions.append(e)

    if len(exceptions) > 0:
        raise exceptions[0]


def _handle_test_output_file(ctx, report_type, test_data, kwds):
    kwd_name = "test_output"
    if report_type != "html":
        kwd_name = "test_output_%s" % report_type

    path = kwds.get(kwd_name, None)
    if path is None:
        message = "No file specified for %s, skipping test output." % kwd_name
        ctx.vlog(message)
        return

    if report_type == "allure":
        file_modication_datatime = kwds.get("file_modication_datatime")
        allure.write_results(path, test_data, file_modication_datatime=file_modication_datatime)
        return

    execution_type = kwds.get("execution_type", "Test")
    try:
        contents = build_report.build_report(test_data, report_type=report_type, execution_type=execution_type)
    except Exception:
        message = f"Problem producing report file {path} for {kwd_name}"
        ctx.vlog(message, exception=True)
        raise

    try:
        with open(path, mode="w", encoding="utf-8") as handle:
            handle.write(unicodify(contents))
    except Exception:
        message = f"Problem writing output file {kwd_name} for {path}"
        ctx.vlog(message, exception=True)
        raise


def _handle_summary(structured_data, **kwds):
    summary_dict = get_dict_value("summary", structured_data)
    num_tests = get_dict_value("num_tests", summary_dict)
    num_failures = get_dict_value("num_failures", summary_dict)
    num_errors = get_dict_value("num_errors", summary_dict)
    num_problems = num_failures + num_errors

    summary_exit_code = EXIT_CODE_OK
    if num_problems > 0:
        summary_exit_code = EXIT_CODE_GENERIC_FAILURE
    elif num_tests == 0:
        summary_exit_code = EXIT_CODE_NO_SUCH_TARGET

    summary_style = kwds.get("summary")
    if kwds.get("test_data_updated"):
        info(TEST_DATA_UPDATED_MESSAGE)
    if summary_style != "none":
        if num_tests == 0:
            warn(NO_TESTS_MESSAGE)
        elif num_problems == 0:
            if kwds.get("update_test_data") and not kwds.get("test_data_updated"):
                info(TEST_DATA_NOT_UPDATED_MESSAGE % num_tests)
            else:
                info(ALL_TESTS_PASSED_MESSAGE % num_tests)
        elif num_problems:
            html_report_file = kwds.get("test_output")
            message_args = (num_problems, num_tests, html_report_file)
            message = PROBLEM_COUNT_MESSAGE % message_args
            warn(message)

        _summarize_tests_full(structured_data, **kwds)

    return summary_exit_code


def _summarize_tests_full(structured_data, **kwds):
    tests = get_dict_value("tests", structured_data)
    for test_case_data in tests:
        _summarize_test_case(test_case_data, **kwds)


def passed(xunit_testcase_el):
    did_pass = True
    for child_el in list(xunit_testcase_el):
        if child_el.tag in ["failure", "error"]:
            did_pass = False
    return did_pass


def _summarize_test_case(structured_data, **kwds):
    summary_style = kwds.get("summary")
    test_id = test_structures.case_id(raw_id=get_dict_value("id", structured_data))
    status = get_dict_value("status", get_dict_value("data", structured_data))
    if status != "success":
        state = click.style("failed", bold=True, fg="red")
    else:
        state = click.style("passed", bold=True, fg="green")
    click.echo(test_id.label + ": " + state)
    if summary_style != "minimal":
        _print_command_line(structured_data, test_id)


def _print_command_line(test, test_id):
    execution_problem = test.get("execution_problem", None)
    if execution_problem:
        click.echo("| command: *could not execute job, no command generated* ")
        return

    job = None
    try:
        job = test["job"]
    except (KeyError, IndexError):
        click.echo("| command: *failed to find job for test object [%s]" % test)
        return
    try:
        command = job["command_line"]
    except (KeyError, IndexError):
        click.echo("| command: *failed to find command_line for job object [%s]" % job)
        return

    click.echo("| command: %s" % command)


__all__ = (
    "handle_reports",
    "handle_reports_and_summary",
)
