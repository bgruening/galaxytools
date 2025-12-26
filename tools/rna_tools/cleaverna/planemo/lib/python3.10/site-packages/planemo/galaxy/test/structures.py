"""Utilities for reasoning about Galaxy test results."""

import os
import shlex
from typing import NamedTuple
from xml.etree import ElementTree as ET

from planemo.io import error
from planemo.test.results import StructuredData as BaseStructuredData

NO_STRUCTURED_FILE = (
    "Warning: Problem with target Galaxy, it did not "
    "produce a structured test results file [%s] - summary "
    "information and planemo reports will be incorrect."
)


class GalaxyTestCommand:
    """Abstraction around building a ``run_tests.sh`` command for Galaxy tests."""

    def __init__(
        self,
        html_report_file,
        xunit_report_file,
        structured_report_file,
        failed=False,
        installed=False,
    ):
        self.html_report_file = html_report_file
        self.xunit_report_file = xunit_report_file
        self.structured_report_file = structured_report_file
        self.failed = failed
        self.installed = installed

    def build(self):
        xunit_report_file = self.xunit_report_file
        sd_report_file = self.structured_report_file
        cmd = "./run_tests.sh $COMMON_STARTUP_ARGS --report_file %s" % shlex.quote(self.html_report_file)
        if xunit_report_file:
            cmd += " --xunit_report_file %s" % shlex.quote(xunit_report_file)
        if sd_report_file:
            cmd += " --structured_data_report_file %s" % shlex.quote(sd_report_file)
        if self.installed:
            cmd += " -installed"
        elif self.failed:
            sd = StructuredData(sd_report_file)
            tests = " ".join(sd.failed_ids)
            cmd += " %s" % tests
        else:
            cmd += " functional.test_toolbox"
        return cmd


class StructuredData(BaseStructuredData):
    """Abstraction around Galaxy's structured test data output."""

    def __init__(self, json_path):
        if not json_path or not os.path.exists(json_path):
            error(NO_STRUCTURED_FILE % json_path)
        super().__init__(json_path)

    def merge_xunit(self, xunit_root):
        self.has_details = True
        xunit_attrib = xunit_root.attrib
        num_tests = int(xunit_attrib.get("tests", 0))
        num_failures = int(xunit_attrib.get("failures", 0))
        num_errors = int(xunit_attrib.get("errors", 0))
        num_skips = int(xunit_attrib.get("skips", 0))
        summary = dict(
            num_tests=num_tests,
            num_failures=num_failures,
            num_errors=num_errors,
            num_skips=num_skips,
        )

        self.structured_data["summary"] = summary

        for testcase_el in xunit_t_elements_from_root(xunit_root):
            test = case_id(testcase_el)
            test_data = self.structured_data_by_id.get(test.id)
            if not test_data:
                continue
            problem_el = None
            for problem_type in ["skip", "failure", "error"]:
                problem_el = testcase_el.find(problem_type)
                if problem_el is not None:
                    break
            if problem_el is not None:
                status = problem_el.tag
                test_data["problem_type"] = problem_el.attrib["type"]
                test_data["problem_log"] = problem_el.text
            else:
                status = "success"
            test_data["status"] = status


class GalaxyTestResults:
    """Class that combine the test-centric xunit output
    with the Galaxy centric structured data output - and
    abstracts away the difference (someday).
    """

    def __init__(
        self,
        output_json_path,
        output_xml_path,
        output_html_path,
        exit_code,
    ):
        self.output_html_path = output_html_path
        sd = StructuredData(output_json_path)
        self.sd = sd
        self.structured_data = sd.structured_data
        self.structured_data_tests = sd.structured_data_tests
        self.structured_data_by_id = sd.structured_data_by_id

        self.xunit_tree = parse_xunit_report(output_xml_path)
        sd.merge_xunit(self._xunit_root)

        self.sd.set_exit_code(exit_code)
        self.sd.read_summary()
        self.sd.update()

    @property
    def exit_code(self):
        return self.sd.exit_code

    @property
    def has_details(self):
        return self.sd.has_details

    @property
    def num_tests(self):
        return self.sd.num_tests

    @property
    def num_problems(self):
        return self.sd.num_problems

    @property
    def _xunit_root(self):
        return self.xunit_tree.getroot()

    @property
    def all_tests_passed(self):
        return self.sd.num_problems == 0

    @property
    def xunit_testcase_elements(self):
        return xunit_t_elements_from_root(self._xunit_root)


def xunit_t_elements_from_root(xunit_root):
    yield from find_cases(xunit_root)


def parse_xunit_report(xunit_report_path):
    return ET.parse(xunit_report_path)


def find_cases(xunit_root):
    return xunit_root.findall("testcase")


def case_id(testcase_el=None, raw_id=None):
    if raw_id is None:
        assert testcase_el is not None
        name_raw = testcase_el.attrib["name"]
        if "TestForTool_" in name_raw:
            raw_id = name_raw
        else:
            class_name = testcase_el.attrib["classname"]
            raw_id = f"{class_name}.{name_raw}"

    name = None
    num = None
    if "TestForTool_" in raw_id:
        tool_and_num = raw_id.split("TestForTool_")[-1]
        if ".test_tool_" in tool_and_num:
            name, num_str = tool_and_num.split(".test_tool_", 1)
            num = _parse_num(num_str)
            # Tempted to but something human friendly in here like
            # num + 1 - but then it doesn't match HTML report.
        else:
            name = tool_and_num
    else:
        name = raw_id

    return TestId(name, num, raw_id)


def _parse_num(num_str):
    try:
        num = int(num_str)
    except ValueError:
        num = None
    return num


class TestId(NamedTuple):
    name: str
    num: int
    id: str

    @property
    def label(self):
        if self.num is not None:
            return f"{self.name}[{self.num}]"
        else:
            return self.id
