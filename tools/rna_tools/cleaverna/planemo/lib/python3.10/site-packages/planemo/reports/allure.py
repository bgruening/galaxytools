import json
import re
import uuid

from allure_commons import plugin_manager
from allure_commons.lifecycle import AllureLifecycle
from allure_commons.logger import AllureFileLogger
from allure_commons.model2 import (
    Label,
    Link,
    Status,
    StatusDetails,
)
from allure_commons.types import (
    AttachmentType,
    LabelType,
    LinkType,
)
from allure_commons.utils import (
    md5,
    platform_label,
    uuid4,
)
from dateutil import parser
from galaxy.util import safe_makedirs

JSON_INDENT = 2
WORKFLOW_INDEX_MATCH = re.compile(r"(.*)\_([\d+])")


class AllureListener:
    def __init__(self, lifecycle):
        self.lifecycle = lifecycle


class AllureWriter:
    def __init__(self, results_path):
        safe_makedirs(results_path)
        self.lifecycle = AllureLifecycle()
        self.logger = AllureFileLogger(results_path)
        self.listener = AllureListener(self.lifecycle)

    def process(self, structured_data, file_modication_datetime=None):
        plugin_manager.register(self.listener)
        plugin_manager.register(self.logger)

        for test_case in structured_data["tests"]:
            self.process_test_case(test_case, file_modication_datetime=file_modication_datetime)

        plugin_manager.unregister(plugin=self.listener)
        plugin_manager.unregister(plugin=self.logger)

    def process_test_case(self, test_case, file_modication_datetime=None):
        with self.lifecycle.schedule_test_case() as test_result:
            test_index = test_case["id"]
            test_data = test_case.get("data") or {}
            job = test_data.get("job") or {}
            test_result.name = test_index
            self._record_start_stop(test_result, file_modication_datetime, job, test_data)

            test_result.fullName = test_index
            test_result.testCaseId = md5(test_index)
            # Maybe an option to swap this - seems like it is causing all runs to seem like
            # retries instead of history proper?
            # test_result.historyId = md5(test_index)
            test_result.historyId = md5(str(uuid.uuid4()))
            tool_id = self._record_suite_labels(test_result, test_index, test_data, job)

            self._attach_data(
                "test_data", json.dumps(test_data, indent=JSON_INDENT), attachment_type=AttachmentType.JSON
            )
            for key in ["stderr", "stdout", "command_line", "external_id", "job_messages"]:
                val = job.get(key)
                if not val:
                    continue
                if isinstance(val, list):
                    attachment_type = AttachmentType.JSON
                    # job messages
                    val = json.dumps(val, indent=JSON_INDENT)
                else:
                    if not val.strip():
                        continue
                    attachment_type = AttachmentType.TEXT
                self._attach_data(key, val, attachment_type=attachment_type)

            problem_message = None
            for key in ["execution_problem", "output_problems"]:
                val = test_data.get(key)
                if not val:
                    continue
                if isinstance(val, list) and val:
                    # remove duplicated messages...
                    val = list(set(val))

                    attachment_type = AttachmentType.HTML
                    as_html_list = "<ul>"
                    as_html_list += "\n".join([f"<li><pre>{v}</pre></li>" for v in val])
                    as_html_list += "</ul>"
                    problem_message = val[0]
                    val = as_html_list
                else:
                    if not val.strip():
                        continue
                    attachment_type = AttachmentType.TEXT
                    problem_message = val
                self._attach_data(key, val, attachment_type=attachment_type)

            if problem_message is None and "job_messages" in job:
                job_messages = job.get("job_messages")
                if job_messages:
                    problem_message = str(job_messages)

            test_result.labels.append(Label(name=LabelType.FRAMEWORK, value="planemo"))
            test_result.labels.append(Label(name=LabelType.LANGUAGE, value=platform_label()))

            self._record_tool_link(test_result, tool_id)
            self._record_status(test_result, test_data)
            if test_result.status in [Status.BROKEN, Status.FAILED]:
                test_result.statusDetails = StatusDetails(message=(problem_message or "Unknown problem"), trace=None)

        self.lifecycle.write_test_case()

    def _record_start_stop(self, test_result, file_modication_datetime, job, test_data):
        start_datetime = file_modication_datetime
        end_datetime = file_modication_datetime
        if "start_datetime" in test_data:
            start_datetime = parser.parse(test_data["start_datetime"])
        elif "create_time" in job:
            start_datetime = parser.parse(job["create_time"])

        if "end_datetime" in test_data:
            end_datetime = parser.parse(test_data["end_datetime"])
        if "update_time" in job:
            end_datetime = parser.parse(job["update_time"])

        if start_datetime is not None:
            test_result.start = int(round(start_datetime.timestamp() * 1000))
        if end_datetime is not None:
            test_result.stop = int(round(end_datetime.timestamp() * 1000))

    def _record_suite_labels(self, test_result, test_index, test_data, job):
        tool_id = None
        index_match = WORKFLOW_INDEX_MATCH.match(test_index)
        if "tool_id" in test_data:
            tool_id = test_data["tool_id"]
            test_result.labels.append(Label(name=LabelType.PARENT_SUITE, value=tool_id))
        elif "tool_id" in job:
            tool_id = job["tool_id"]
            test_result.labels.append(Label(name=LabelType.PARENT_SUITE, value=tool_id))
        elif index_match:
            test_result.labels.append(Label(name=LabelType.PARENT_SUITE, value=index_match.group(1)))

        if "tool_version" in test_data:
            test_result.labels.append(Label(name=LabelType.SUITE, value=test_data["tool_version"]))
        elif "tool_version" in job:
            test_result.labels.append(Label(name=LabelType.SUITE, value=job["tool_version"]))

        if "test_index" in test_data:
            test_result.labels.append(Label(name=LabelType.SUB_SUITE, value=str(test_data["test_index"])))
        elif index_match:
            test_result.labels.append(Label(name=LabelType.SUB_SUITE, value=index_match.group(2)))

        return tool_id

    def _record_tool_link(self, test_result, tool_id):
        if tool_id and "repos" in tool_id:
            tool_parts = tool_id.split("/")
            if len(tool_parts) >= 4:
                link = Link(LinkType.LINK, "https://%s" % "/".join(tool_parts[0:4]), "Tool Repository")
                test_result.links.append(link)

    def _record_status(self, test_result, test_data):
        status = test_data.get("status", "error")
        if status == "success":
            test_result.status = Status.PASSED
        elif status == "failure":
            test_result.status = Status.FAILED
        elif status == "skip":
            test_result.status = Status.SKIPPED
        else:
            test_result.status = Status.BROKEN

    def _attach_data(self, key, val, attachment_type=AttachmentType.TEXT):
        self.lifecycle.attach_data(uuid4(), val, name=key, attachment_type=attachment_type, extension=None)


def write_results(results_path, structured_data, **kwds):
    AllureWriter(results_path).process(structured_data)
