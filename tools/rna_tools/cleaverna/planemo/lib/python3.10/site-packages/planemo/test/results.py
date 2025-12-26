"""Describes results.

Is a JSON.
"""

import json
import os

from planemo.io import error


class StructuredData:
    """Abstraction around a simple data structure describing test results."""

    def __init__(self, json_path=None, data=None):
        """Create a :class:`StructuredData` from a JSON file."""

        def data_error():
            error(
                "An invalid JSON for structured test result data - "
                "summary information and planemo reports will be "
                "incorrect."
            )

        self.json_path = json_path
        structured_data = {}
        structured_data_tests = {}
        if json_path and os.path.exists(json_path) and data is None:
            try:
                with open(json_path) as output_json_f:
                    data = json.load(output_json_f)
            except Exception:
                data_error()

        try:
            structured_data = data
            structured_data_tests = structured_data["tests"]
        except Exception:
            data_error()

        self.structured_data = structured_data
        self.structured_data_tests = structured_data_tests
        structured_data_by_id = {}
        for test in self.structured_data_tests:
            structured_data_by_id[test["id"]] = test["data"]
        self.structured_data_by_id = structured_data_by_id
        self.has_details = "summary" in structured_data
        if self.has_details:
            self.read_summary()

    def update(self):
        """Write out an updated version of this data structure to supplied json path."""
        with open(self.json_path, "w") as out_f:
            json.dump(self.structured_data, out_f)

    def set_exit_code(self, exit_code):
        """Set the exit_code for the this test."""
        self.structured_data["exit_code"] = exit_code

    def calculate_summary_data_if_needed(self):
        if "summary" not in self.structured_data:
            self.calculate_summary_data()

    def calculate_summary_data(self):
        """Use full details on individual test data to update structured data with summary info."""
        num_tests = 0
        num_failures = 0
        num_skips = 0
        num_errors = 0

        for test in self.structured_data_tests:
            test_data = get_dict_value("data", test)
            status = get_dict_value("status", test_data)
            num_tests += 1
            if status == "skip":
                num_skips += 1
            elif status == "failure":
                num_failures += 1
            elif status == "error":
                num_errors += 1
            elif status != "success":
                raise Exception("Unknown test status encountered [%s]" % status)

        summary = {}
        summary["num_tests"] = num_tests
        summary["num_failures"] = num_failures
        summary["num_skips"] = num_skips
        summary["num_errors"] = num_errors
        self.structured_data["summary"] = summary
        self.read_summary()

    def read_summary(self):
        """Read summary data into properties on this class."""
        summary = self.structured_data["summary"]
        num_tests = summary["num_tests"]
        num_failures = summary["num_failures"]
        num_skips = summary["num_skips"]
        num_errors = summary["num_errors"]

        self.num_tests = num_tests
        self.num_problems = num_skips + num_errors + num_failures

        self.exit_code = self.structured_data.get("exit_code", None)

    @property
    def failed_ids(self):
        """Find set of IDs for failed tests."""
        ids = set()
        for test_data in self.structured_data_tests:
            if test_data["data"]["status"] == "success":
                continue
            test_case = test_data["id"].replace(".test_toolbox.", ".test_toolbox:")
            ids.add(test_case)
        return ids


def get_dict_value(key, data):
    """Return data[key] with improved KeyError."""
    try:
        return data[key]
    except (KeyError, TypeError):
        raise KeyError(f"No key [{key}] in [{data}]")


__all__ = (
    "StructuredData",
    "get_dict_value",
)
