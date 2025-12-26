import base64

from galaxy.util import strip_control_characters
from galaxy.util.resources import resource_string
from jinja2 import (
    Environment,
    PackageLoader,
)

TITLE = "Results (powered by Planemo)"

cancel_fragment = "Invocation scheduling cancelled because"
fail_fragment = "Invocation scheduling failed because"


def render_message_to_string(invocation_message):  # noqa: C901
    # ChatGPT did a reasonable job of translating this from https://github.com/galaxyproject/galaxy/blob/d92bbb144ffcda7e17368cf43dd25c8a9a3a7dd6/client/src/components/WorkflowInvocationState/InvocationMessage.vue#L93-L172
    reason = invocation_message["reason"]
    if reason == "user_request":
        return f"{cancel_fragment} user requested cancellation."
    elif reason == "history_deleted":
        return f"{cancel_fragment} the history of the invocation was deleted."
    elif reason == "cancelled_on_review":
        return f"{cancel_fragment} the invocation was paused at step {invocation_message['workflow_step_id'] + 1} and not approved."
    elif reason == "collection_failed":
        return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} requires a dataset collection created by step {invocation_message['dependent_workflow_step_id'] + 1}, but dataset collection entered a failed state."
    elif reason == "dataset_failed":
        if invocation_message.get("dependent_workflow_step_id") is not None:
            return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} requires a dataset created by step {invocation_message['dependent_workflow_step_id'] + 1}, but dataset entered a failed state."
        else:
            return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} requires a dataset, but dataset entered a failed state."
    elif reason == "job_failed":
        return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} depends on job(s) created in step {invocation_message['dependent_workflow_step_id'] + 1}, but a job for that step failed."
    elif reason == "output_not_found":
        return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} depends on output '{invocation_message['output_name']}' of step {invocation_message['dependent_workflow_step_id'] + 1}, but this step did not produce an output of that name."
    elif reason == "expression_evaluation_failed":
        return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} contains an expression that could not be evaluated."
    elif reason == "when_not_boolean":
        return f"{fail_fragment} step {invocation_message['workflow_step_id'] + 1} is a conditional step and the result of the when expression is not a boolean type."
    elif reason == "unexpected_failure":
        at_step = ""
        if invocation_message.get("workflow_step_id") is not None:
            at_step = f" at step {invocation_message['workflow_step_id'] + 1}"
        if "details" in invocation_message and invocation_message["details"]:
            return f"{fail_fragment} an unexpected failure occurred{at_step}: '{invocation_message['details']}'"
        return f"{fail_fragment} an unexpected failure occurred{at_step}."
    elif reason == "workflow_output_not_found":
        return f"Defined workflow output '{invocation_message['output_name']}' was not found in step {invocation_message['workflow_step_id'] + 1}."
    else:
        return reason


def build_report(structured_data, report_type="html", execution_type="Test", **kwds):
    """Use report_{report_type}.tpl to build page for report."""
    environment = dict(
        title=TITLE,
        raw_data=structured_data,
    )

    __fix_test_ids(environment)
    environment = __inject_summary(environment)
    environment["execution_type"] = execution_type

    if report_type == "html":
        # The HTML report format needs a lot of extra, custom data.
        # IMO, this seems to suggest it should be embedded.
        environment["title"] = None
        markdown = template_data(environment, "report_markdown.tpl")
        environment["title"] = " ".join((environment["execution_type"], TITLE))
        environment["raw_data"] = base64.b64encode(markdown.encode("utf-8")).decode("utf-8")
        environment.update(
            {
                "custom_style": __style("custom.css"),
                "custom_script": __script("custom"),
                "bootstrap_style": __style("bootstrap.min.css"),
                "jquery_script": __script("jquery.min"),
                "bootstrap_script": __script("bootstrap.min"),
                "markdown_it_script": __script("markdown-it.min"),
            }
        )

    return template_data(environment, "report_%s.tpl" % report_type)


def template_data(environment, template_name, **kwds):
    """Build an arbitrary templated page."""
    env_kwargs = {}
    if template_name == "report_markdown.tpl":
        env_kwargs["keep_trailing_newline"] = True
        env_kwargs["trim_blocks"] = True
    env = Environment(loader=PackageLoader("planemo", "reports"), **env_kwargs)
    env.filters["strip_control_characters"] = lambda x: strip_control_characters(x) if x else x
    env.globals["render_message_to_string"] = render_message_to_string
    template = env.get_template(template_name)
    return template.render(**environment)


def __fix_test_ids(environment):
    for test in environment["raw_data"]["tests"]:
        test_data = test.get("data")
        if test_data and test_data.get("tool_id"):
            test["id"] = f"{test_data['tool_id']} (Test #{test_data['test_index'] + 1})"


def __inject_summary(environment):
    total = 0
    errors = 0
    failures = 0
    skips = 0
    for execution in environment["raw_data"]["tests"]:
        total += 1
        test_data = execution.get("data")
        if test_data:
            status = test_data.get("status")
            if status == "error":
                errors += 1
            elif status == "failure":
                failures += 1
            elif status == "skipped":
                skips += 1
    environment["raw_data"]["results"] = {
        "total": total,
        "errors": errors,
        "failures": failures,
        "skips": skips,
    }
    if "suitename" not in environment:
        environment["raw_data"]["suitename"] = TITLE
    return environment


def __style(filename):
    resource = resource_string("planemo.reports", filename)
    return "<style>%s</style>" % resource


def __script(short_name):
    resource = resource_string("planemo.reports", "%s.js" % short_name)
    return "<script>%s</script>" % resource
