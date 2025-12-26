import random
from io import StringIO
from typing import (
    Dict,
    List,
    Optional,
    Set,
)

from rich.console import Group
from rich.live import Live
from rich.panel import Panel
from rich.progress import (
    BarColumn,
    Progress,
    TaskID,
    TaskProgressColumn,
    TextColumn,
)
from typing_extensions import TypedDict

from .api import (
    invocation_state_terminal,
    InvocationApi,
    JOB_ERROR_STATES,
)
from .progress_display import DisplayConfiguration


# Types for various invocation responses
class InvocationStep(TypedDict, total=False):
    id: str
    state: Optional[str]
    subworkflow_invocation_id: Optional[str]


class Invocation(TypedDict, total=False):
    id: str
    state: str
    steps: List[InvocationStep]


class InvocationJobsSummary(TypedDict, total=False):
    states: Dict[str, int]


class WorkflowProgress(Progress):
    _jobs_task: TaskID
    _steps_task: TaskID
    _subworkflows_task: Optional[TaskID] = None

    def __init__(self, display: DisplayConfiguration):
        self.invocation_state: str = "new"
        self.step_count: Optional[int] = None
        self.job_count: Optional[int] = 0
        self.jobs_completed: Optional[int] = None
        self.step_states: Dict[str, int] = {}
        self.num_ok: int = 0
        self.num_new: int = 0
        self.num_queued: int = 0
        self.num_running: int = 0
        self.num_errors: int = 0
        self.num_paused: int = 0
        self.printed_job_errors: Set[str] = set()  # Track printed job errors by job ID

        self.num_subworkflows: int = 0
        self.num_subworkflows_complete: int = 0
        self.display = display
        bar_column = BarColumn(
            style=self.display.style_bar_back,
            finished_style=self.display.style_bar_finished,
            complete_style=self.display.style_bar_complete,
        )
        self.jobs_color: str = self.display.style_initializing
        self.steps_color: str = self.display.style_initializing
        self.subworkflows_color: str = self.display.style_initializing
        super().__init__(
            TextColumn("[progress.description]{task.description}"),
            TextColumn(display.divider),
            bar_column,
            TextColumn(display.divider),
            TaskProgressColumn(f"[{self.display.style_percent}]{{task.percentage:>3.0f}}%"),
            TextColumn(display.divider),
            TextColumn(text_format="{task.fields[status]}"),
        )
        self.add_bars()

    @property
    def invocation_scheduling_terminal(self):
        # This is a tricky issue, because we might be waiting to schedule a new step whose inputs are paused.
        # Unlikely that this will ever be unpaused, so we also consider the "ready" invocation state and the presence of paused jobs as terminal.
        # We might want to tweak the state in Galaxy if we pause jobs because of input errors, so that this hack won't be rquired.
        return (
            invocation_state_terminal(self.invocation_state)
            or self.invocation_state == "ready"
            and self.num_paused
            and self.num_errors
        )

    @property
    def jobs_terminal(self):
        return self.job_count is not None and self.job_count == self.jobs_terminal_count

    @property
    def terminal(self):
        return self.invocation_scheduling_terminal and self.jobs_terminal

    def handle_subworkflow_counts(self, num: int, num_complete: int):
        previous_count = self.num_subworkflows
        self.num_subworkflows = num
        self.num_subworkflows_complete = num_complete
        if previous_count < 2 and num >= 2:
            self._subworkflows_task = self.add_task(
                f"[{self.subworkflows_color}]{self.display.label_progress_subworkflows}", status=""
            )

        if num >= 2:
            self.subworkflows_color = self.display.style_ok
            subworkflows_status = f"{self.num_subworkflows_complete}/{self.num_subworkflows} terminal"
            self.update(
                self._subworkflows_task,
                total=self.num_subworkflows,
                completed=self.num_subworkflows_complete,
                description=f"[{self.subworkflows_color}]{self.display.label_progress_subworkflows}",
                status=subworkflows_status,
            )

    def handle_invocation(self, invocation: Invocation, job_state_summary: InvocationJobsSummary):
        self.invocation_state = invocation.get("state") or "new"
        self.step_count = len(invocation.get("steps") or []) or None
        self.step_states = step_states(invocation)

        steps_completed = None

        steps_status = ""
        if self.step_count is None:
            steps_status = "Loading steps."
            self.steps_color = self.display.style_initializing
        elif self.invocation_state == "cancelled":
            steps_status = "Invocation cancelled"
            self.steps_color = self.display.style_error
        elif self.invocation_state == "failed":
            steps_status = "Invocation failed"
            self.steps_color = self.display.style_error
        else:
            num_scheduled = self.step_states.get("scheduled") or 0
            if num_scheduled > 0:
                self.steps_color = self.display.style_ok
            else:
                self.steps_color = self.display.style_initializing
            steps_completed = num_scheduled
            steps_status = f"{num_scheduled}/{self.step_count} scheduled"

        jobs_status = ""
        self.job_count = job_count(job_state_summary)
        self.num_new = count_states(job_state_summary, ["new"])
        self.num_queued = count_states(job_state_summary, ["queued", "waiting"])
        self.num_running = count_states(job_state_summary, ["running"])
        self.num_errors = error_count(job_state_summary)
        self.num_ok = ok_count(job_state_summary)
        self.jobs_completed = self.num_ok + self.num_errors
        self.num_paused = count_states(job_state_summary, ["paused"])
        self.jobs_terminal_count = self.jobs_completed + self.num_paused
        jobs_total: Optional[int] = self.job_count
        if self.num_errors > 0:
            self.jobs_color = self.display.style_error
        elif self.job_count > 0:
            self.jobs_color = self.display.style_ok
        else:
            self.jobs_color = self.display.style_initializing
            self.jobs_completed = None
            jobs_total = None
        if self.job_count > 0:
            jobs_status = f"{self.jobs_completed}/{self.job_count} terminal"
        self.update(
            self._steps_task,
            total=self.step_count,
            completed=steps_completed,
            description=f"[{self.steps_color}]{self.display.label_progress_steps}",
            status=steps_status,
        )
        self.update(
            self._jobs_task,
            total=jobs_total,
            completed=self.jobs_completed,
            description=f"[{self.jobs_color}]{self.display.label_progress_jobs}",
            status=jobs_status,
        )

    def _job_states_console_line(self):
        output = StringIO()
        if self.num_ok > 0:
            output.write(f"{self.display.icon_state_ok} {self.num_ok} {self.display.divider} ")
        if self.num_errors > 0:
            output.write(f"{self.display.icon_state_errors} {self.num_errors} {self.display.divider} ")
        if self.num_new > 0:
            output.write(f"{self.display.icon_state_new} {self.num_new} {self.display.divider} ")
        if self.num_queued > 0:
            output.write(f"{self.display.icon_state_queued} {self.num_queued} {self.display.divider} ")
        if self.num_running > 0:
            output.write(f"{self.display.icon_state_running} {self.num_running} {self.display.divider} ")
        if self.num_paused > 0:
            output.write(f"{self.display.icon_state_paused} {self.num_paused} {self.display.divider} ")

        result = output.getvalue().rstrip(" {self.display.divider} ")
        output.close()
        # Is there an actual way to reset this? The undefined style seems to work but is a hack.
        return f"[{self.jobs_color}]{self.display.label_job_states_prefix} [reset]{self.display.divider} {result}"

    def add_bars(self):
        self._steps_task = self.add_task(f"[{self.steps_color}]{self.display.label_progress_steps}", status="")
        self._jobs_task = self.add_task(f"[{self.jobs_color}]{self.display.label_progress_jobs}", status="")

    def print_job_errors_once(
        self,
        ctx,
        invocation_api: InvocationApi,
        invocation_id: str,
        workflow_progress_display: "WorkflowProgressDisplay",
    ):
        """Print job errors only if they haven't been printed before, tracking by job ID."""

        # Early exit if no errors detected in current state
        if self.num_errors == 0:
            return

        # Get all failed job details and filter for new ones
        try:
            new_error_lines = []

            # Use the efficient jobs API to get all failed jobs for this invocation in one call
            try:
                # Get all jobs in error state for this invocation
                failed_jobs = invocation_api.get_invocation_jobs(invocation_id, state="error")
                ctx.vlog(f"Found {len(failed_jobs)} failed jobs for invocation {invocation_id}")

                # Process each failed job that we haven't seen before
                for job in failed_jobs:
                    job_id = job.get("id")
                    if job_id and job_id not in self.printed_job_errors:
                        # Mark this job as printed
                        self.printed_job_errors.add(job_id)
                        ctx.vlog(f"Processing failed job {job_id}")

                        # Get detailed job information for error details
                        try:
                            job_details = invocation_api.get_job(job_id, full_details=True)
                            job_errors = self._format_job_error_details(job_id, job_details)
                            new_error_lines.extend(job_errors)
                        except Exception as e:
                            ctx.vlog("Failed to get details for job {}: {}".format(job_id, e))

            except Exception as e:
                ctx.vlog("Failed to get failed jobs: {}".format(e))
                return

            # Print any new errors found
            if new_error_lines:
                workflow_progress_display.console.print("\n".join(new_error_lines))

        except Exception as e:
            error_msg = "❌ Failed to collect failed job details: {}".format(e)
            workflow_progress_display.console.print(error_msg)

    def _format_job_error_details(self, job_id, job_details):
        """Format error details for a single failed job."""
        error_lines = []
        exit_code = job_details.get("exit_code")
        stderr = job_details.get("stderr", "").strip()
        stdout = job_details.get("stdout", "").strip()
        command_line = job_details.get("command_line")
        tool_id = job_details.get("tool_id")

        error_lines.append("❌ Failed job {}:".format(job_id))
        if tool_id:
            error_lines.append("   Tool ID: {}".format(tool_id))
        if exit_code is not None:
            error_lines.append("   Exit code: {}".format(exit_code))
        if command_line:
            error_lines.append("   Command line: {}".format(command_line))
        if stdout:
            error_lines.append("   Stdout:")
            for line in stdout.split("\n"):
                if line.strip():
                    error_lines.append("     {}".format(line))
        if stderr:
            error_lines.append("   Stderr:")
            for line in stderr.split("\n"):
                if line.strip():
                    error_lines.append("     {}".format(line))
        error_lines.append("")  # Empty line for spacing
        return error_lines


# converted from Galaxy TypeScript (see util.ts next to WorkflowInvocationState.vue)
def count_states(job_summary: Optional[InvocationJobsSummary], query_states: list[str]) -> int:
    count = 0
    states = job_summary.get("states") if job_summary else None
    if states:
        for state in query_states:
            count += states.get(state, 0)
    return count


def job_count(job_summary: Optional[InvocationJobsSummary]) -> int:
    states = job_summary.get("states") if job_summary else None
    count = 0
    if states:
        for state_count in states.values():
            if state_count:
                count += state_count
    return count


def step_states(invocation: Invocation):
    step_states = {}
    steps = invocation.get("steps") or []
    for step in steps:
        if not step:
            continue
        step_state = step.get("state") or "unknown"
        if step_state not in step_states:
            step_states[step_state] = 0
        step_states[step_state] += 1

    return step_states


def ok_count(job_summary: InvocationJobsSummary) -> int:
    return count_states(job_summary, ["ok", "skipped"])


def error_count(job_summary: InvocationJobsSummary) -> int:
    return count_states(job_summary, JOB_ERROR_STATES)


def running_count(job_summary: InvocationJobsSummary) -> int:
    return count_states(job_summary, ["running"])


class WorkflowProgressDisplay(Live):
    def __init__(
        self,
        invocation_id: str,
        display_configuration: Optional[DisplayConfiguration] = None,
        galaxy_url: Optional[str] = None,
    ):
        self.subworkflow_invocation_ids_seen: Set[str] = set()
        self.subworkflow_invocation_ids_completed: Set[str] = set()
        self.subworkflow_invocation_id: Optional[str] = None
        self.new_steps: List[str] = []
        self.invocation_id = invocation_id
        display = display_configuration or DisplayConfiguration()
        self.galaxy_url = galaxy_url
        self.display = display
        self.workflow_progress = WorkflowProgress(display)
        self.subworkflow_progress = WorkflowProgress(display)
        super().__init__(self._panel(), auto_refresh=False)

    def _register_subworkflow_invocation_ids_from(self, invocation: Invocation):
        subworkflow_invocation_ids: List[str] = []
        steps = invocation.get("steps") or []
        new_steps: List[str] = []
        for step in steps:
            if step["state"] == "new":
                new_steps.append(step["id"])
            subworkflow_invocation_id = step.get("subworkflow_invocation_id")
            if subworkflow_invocation_id:
                subworkflow_invocation_ids.append(subworkflow_invocation_id)
        self.new_steps = new_steps
        self._register_subworkflow_invocation_ids(subworkflow_invocation_ids)

    def _register_subworkflow_invocation_ids(self, ids: List[str]):
        for invocation_id in ids:
            self.subworkflow_invocation_ids_seen.add(invocation_id)

    def _complete_subworkflow(self, id: str):
        self.subworkflow_invocation_ids_completed.add(id)

    def an_incomplete_subworkflow_id(self):
        return random.choice(tuple(self.subworkflow_invocation_ids_seen - self.subworkflow_invocation_ids_completed))

    def all_subworkflows_complete(self):
        if self.new_steps:
            # These don't have subworkflow invocation ids yet, we can't know if they're all complete
            return False
        return len(self.subworkflow_invocation_ids_seen) == len(self.subworkflow_invocation_ids_completed)

    def get_invocation_ui_link(self):
        if self.galaxy_url:
            return f"{self.galaxy_url}/workflows/invocations/{self.invocation_id}"
        else:
            return None

    def _panel(self):
        def job_states(workflow_progress):
            if self.display.include_job_state_breakdown:
                return workflow_progress._job_states_console_line()
            else:
                return None

        title = f"[{self.display.style_header}]{self.display.label_header_prefix}<[link={self.get_invocation_ui_link()}]{self.invocation_id}[/link]>"
        subworkflow_title = None
        if self.subworkflow_invocation_id:
            subworkflow_title = f"[{self.display.style_subworkflow_header}]{self.display.label_subworkflow_header_prefix}<{self.subworkflow_invocation_id}>"

        if not self.subworkflow_invocation_id or not self.display.include_nested_subworkflows:
            renderable = as_group(
                self.workflow_progress,
                job_states(self.workflow_progress),
            )
        elif not self.display.subworkflows_as_panel:
            renderable = as_group(
                self.workflow_progress,
                job_states(self.workflow_progress),
                subworkflow_title,
                self.subworkflow_progress,
                job_states(self.subworkflow_progress),
            )
        else:
            renderable = as_group(
                self.workflow_progress,
                job_states(self.workflow_progress),
                Panel(
                    as_group(self.subworkflow_progress, job_states(self.subworkflow_progress)),
                    title=subworkflow_title,
                ),
            )
        return Panel(renderable, title=title, expand=True)

    def _update_panel(self):
        self.update(self._panel())
        self.refresh()

    def handle_invocation(self, invocation: Invocation, job_state_summary: InvocationJobsSummary):
        self.workflow_progress.handle_invocation(invocation, job_state_summary)
        self._register_subworkflow_invocation_ids_from(invocation)
        self._update_panel()

    def handle_subworkflow_invocation(self, invocation: Invocation, job_state_summary: InvocationJobsSummary):
        self.subworkflow_invocation_id = invocation["id"]
        self.subworkflow_progress.handle_invocation(invocation, job_state_summary)
        self._register_subworkflow_invocation_ids_from(invocation)
        if self.subworkflow_progress.terminal:
            self._complete_subworkflow(invocation["id"])
        self.workflow_progress.handle_subworkflow_counts(
            len(self.subworkflow_invocation_ids_seen),
            len(self.subworkflow_invocation_ids_completed),
        )
        self._update_panel()


def as_group(*renderables):
    return Group(*(r for r in renderables if r))
