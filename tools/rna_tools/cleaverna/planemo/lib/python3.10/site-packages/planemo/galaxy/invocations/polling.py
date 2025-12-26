import time
from typing import (
    List,
    Optional,
    Protocol,
)

from .api import (
    invocation_state_terminal,
    InvocationApi,
    InvocationJobsSummary,
    JOB_ERROR_STATES,
    NON_TERMINAL_JOB_STATES,
)
from .progress import WorkflowProgressDisplay


class PollingTracker(Protocol):
    def sleep(self) -> None: ...


class PollingTrackerImpl(PollingTracker):
    def __init__(self, polling_backoff: int, timeout=None):
        self.polling_backoff = polling_backoff
        self.timeout = timeout
        self.delta = 0.25
        self.total_wait_time = 0

    def sleep(self):
        if self.timeout is not None and self.total_wait_time > self.timeout:
            message = "Timed out while polling Galaxy."
            raise Exception(message)
        self.total_wait_time += self.delta
        time.sleep(self.delta)
        self.delta += self.polling_backoff


def _summarize_invocation(invocation_api: InvocationApi, invocation_id: str):
    invocation = invocation_api.get_invocation(invocation_id)
    assert invocation
    invocation_jobs = invocation_api.get_invocation_summary(invocation_id)
    return invocation, invocation_jobs


def _poll_main_workflow(
    ctx,
    invocation_id: str,
    invocation_api: InvocationApi,
    workflow_progress_display: WorkflowProgressDisplay,
    fail_fast: bool,
):
    if workflow_progress_display.workflow_progress.terminal:
        return None, None, None

    try:
        invocation, invocation_jobs = _summarize_invocation(invocation_api, invocation_id)
        workflow_progress_display.handle_invocation(invocation, invocation_jobs)
        return invocation, invocation_jobs, None
    except Exception as e:
        print(e)
        return None, None, e


def _poll_subworkflow(
    ctx,
    invocation_id: str,
    invocation_api: InvocationApi,
    workflow_progress_display: WorkflowProgressDisplay,
    fail_fast: bool,
):
    if workflow_progress_display.all_subworkflows_complete():
        return None, None, None

    try:
        subworkflow_id = workflow_progress_display.an_incomplete_subworkflow_id()
        invocation, invocation_jobs = _summarize_invocation(invocation_api, subworkflow_id)
        workflow_progress_display.handle_subworkflow_invocation(invocation, invocation_jobs)
        return invocation, invocation_jobs, None
    except Exception as e:
        return None, None, e


def _check_for_errors(
    ctx,
    invocation_id: str,
    exception: Optional[Exception],
    invocation,
    invocation_jobs,
    invocation_api: InvocationApi,
    workflow_progress_display: WorkflowProgressDisplay,
    fail_fast: bool,
):
    error_message = workflow_in_error_message(
        ctx,
        invocation_id,
        exception,
        invocation,
        invocation_jobs,
        invocation_api=invocation_api,
        workflow_progress_display=workflow_progress_display,
        fail_fast=fail_fast,
    )
    if error_message:
        final_state = "new" if not invocation else invocation["state"]
        job_state = summary_job_state(invocation_jobs, fail_fast)
        return final_state, job_state, error_message
    return None


def _is_polling_complete(workflow_progress_display: WorkflowProgressDisplay) -> bool:
    return (
        workflow_progress_display.workflow_progress.terminal and workflow_progress_display.all_subworkflows_complete()
    )


def wait_for_invocation_and_jobs(
    ctx,
    invocation_id: str,
    invocation_api: InvocationApi,
    polling_tracker: PollingTracker,
    workflow_progress_display: WorkflowProgressDisplay,
    fail_fast: bool = False,
):
    ctx.vlog("Waiting for invocation [%s]" % invocation_id)

    last_invocation = None
    last_invocation_jobs = None
    error_message: Optional[str] = None

    while not _is_polling_complete(workflow_progress_display):
        # Poll main workflow
        main_invocation, main_jobs, main_exception = _poll_main_workflow(
            ctx, invocation_id, invocation_api, workflow_progress_display, fail_fast
        )

        if main_invocation:
            last_invocation = main_invocation
            last_invocation_jobs = main_jobs

        error_result = _check_for_errors(
            ctx,
            invocation_id,
            main_exception,
            main_invocation,
            main_jobs,
            invocation_api=invocation_api,
            workflow_progress_display=workflow_progress_display,
            fail_fast=fail_fast,
        )
        if error_result:
            return error_result

        # Poll subworkflow
        sub_invocation, sub_jobs, sub_exception = _poll_subworkflow(
            ctx, invocation_id, invocation_api, workflow_progress_display, fail_fast
        )

        if sub_invocation:
            error_result = _check_for_errors(
                ctx,
                sub_invocation["id"] if sub_invocation else invocation_id,
                sub_exception,
                sub_invocation,
                sub_jobs,
                invocation_api,
                workflow_progress_display,
                fail_fast,
            )
            if error_result:
                return error_result

        if not _is_polling_complete(workflow_progress_display):
            polling_tracker.sleep()

    ctx.vlog(f"The final state of all jobs and subworkflow invocations for invocation [{invocation_id}] is 'ok'")
    job_state = summary_job_state(last_invocation_jobs, fail_fast)
    assert last_invocation

    # Final check for job errors when fail_fast is enabled
    if fail_fast and job_state in JOB_ERROR_STATES and not error_message:
        error_message = workflow_in_error_message(
            ctx,
            invocation_id,
            None,
            last_invocation,
            last_invocation_jobs,
            fail_fast=fail_fast,
            invocation_api=invocation_api,
            workflow_progress_display=workflow_progress_display,
        )

    return last_invocation["state"], job_state, error_message


def workflow_in_error_message(
    ctx,
    invocation_id,
    last_exception,
    last_invocation,
    last_invocation_jobs,
    invocation_api: InvocationApi,
    workflow_progress_display: WorkflowProgressDisplay,
    fail_fast=False,
) -> Optional[str]:
    """Return an error message if workflow is in an error state."""

    invocation_state = "new" if not last_invocation else last_invocation["state"]
    job_state = summary_job_state(last_invocation_jobs, fail_fast)

    error_message = None
    if last_exception:
        ctx.vlog(f"Problem waiting on invocation: {str(last_exception)}")
        error_message = f"Final state of invocation {invocation_id} is [{invocation_state}]"

    if invocation_state_terminal(invocation_state) and invocation_state != "scheduled":
        msg = f"Failed to run workflow, invocation ended in [{invocation_state}] state."
        ctx.vlog(msg)
        error_message = msg if not error_message else f"{error_message}. {msg}"

    # Print job errors when detected, regardless of fail_fast setting
    if job_state in JOB_ERROR_STATES:
        # Print failed job details when we detect job failures, using WorkflowProgress to avoid duplicates
        if invocation_api and workflow_progress_display:
            # Pass the Live display to print errors above the live panel
            workflow_progress_display.workflow_progress.print_job_errors_once(
                ctx, invocation_api, invocation_id, workflow_progress_display=workflow_progress_display
            )

        # Only return error message (which stops execution) when fail_fast is enabled
        if fail_fast:
            msg = f"Failed to run workflow, at least one job is in [{job_state}] state."
            ctx.vlog(msg)
            error_message = msg if not error_message else f"{error_message}. {msg}"

    return error_message


def summary_job_state(job_states_summary: Optional[InvocationJobsSummary], fail_fast: bool = False):
    states = {state for state in (job_states_summary or {"states": {}})["states"]}
    if not fail_fast:
        current_non_terminal_states = NON_TERMINAL_JOB_STATES.intersection(states)
        if current_non_terminal_states:
            # ensure all non-terminal states advance, then return the first failing state, if any.
            return next(iter(current_non_terminal_states))
    if states:
        # We have ensured that that all jobs are terminal, we want to return failed jobs in the summary if there are any.
        for error_state in JOB_ERROR_STATES:
            if error_state in states:
                return error_state
        return next(iter(states))
    else:
        return "ok"


def subworkflow_invocation_ids(invocation_api: InvocationApi, invocation_id: str) -> List[str]:
    invocation = invocation_api.get_invocation(invocation_id)
    subworkflow_invocation_ids = []
    for step in invocation["steps"]:
        subworkflow_invocation_id = step.get("subworkflow_invocation_id")
        if subworkflow_invocation_id:
            subworkflow_invocation_ids.append(subworkflow_invocation_id)
    return subworkflow_invocation_ids
