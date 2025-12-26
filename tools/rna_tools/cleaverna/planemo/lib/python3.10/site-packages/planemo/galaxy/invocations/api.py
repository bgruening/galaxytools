"""API interaction for Galaxy's workflow invocation API.

Gives a mockable surface for testing, type contract consumed by Planemo,
and builtin utilities around bioblend for working around transient request
issues that have been observed in practice.
"""

from typing import (
    Dict,
    List,
    Optional,
    Protocol,
)

from typing_extensions import TypedDict

from planemo.galaxy.api import retry_on_timeouts


class InvocationStep(TypedDict, total=False):
    state: Optional[str]
    subworkflow_invocation_id: Optional[str]


class Invocation(TypedDict, total=False):
    id: str
    state: str
    steps: List[InvocationStep]
    history_id: Optional[str]


class InvocationJobsSummary(TypedDict, total=False):
    states: Dict[str, int]


class Job(TypedDict, total=False):
    id: str
    state: str
    exit_code: Optional[int]
    stderr: Optional[str]
    stdout: Optional[str]
    command_line: Optional[str]
    tool_id: Optional[str]
    history_id: Optional[str]


class InvocationApi(Protocol):
    def get_invocation(self, invocation_id: str) -> Invocation: ...

    def get_invocation_summary(self, invocation_id: str, state: Optional[str] = None) -> InvocationJobsSummary: ...

    def get_invocation_jobs(self, invocation_id: str, state: Optional[str] = None) -> List[Job]: ...

    def get_job(self, job_id: str, full_details: bool = False) -> Job: ...


class BioblendInvocationApi(InvocationApi):
    def __init__(self, ctx, user_gi):
        self._ctx = ctx
        self._user_gi = user_gi

    def get_invocation(self, invocation_id: str) -> Invocation:
        return retry_on_timeouts(self._ctx, self._user_gi, lambda gi: gi.invocations.show_invocation(invocation_id))

    def get_invocation_summary(self, invocation_id: str, state: Optional[str] = None) -> InvocationJobsSummary:
        return self._user_gi.invocations.get_invocation_summary(invocation_id)

    def get_invocation_jobs(self, invocation_id: str, state: Optional[str] = None) -> List[Job]:
        return self._user_gi.jobs.get_jobs(invocation_id=invocation_id, state=state)

    def get_job(self, job_id: str, full_details: bool = False) -> Job:
        return self._user_gi.jobs.show_job(job_id, full_details=full_details)


def invocation_state_terminal(state: str):
    return state in ["scheduled", "cancelled", "failed"]


JOB_ERROR_STATES = ["error", "deleted", "failed", "stopped", "stop", "deleting"]
NON_TERMINAL_JOB_STATES = {"running", "queued", "new", "ready", "resubmitted", "upload", "waiting"}
