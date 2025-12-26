from pydantic import BaseModel

# from rich.style import StyleType - doesn't work with Pydantic, keeping to string styles for now
StyleType = str

# uses a bit more space but has better visual separation between subworkflows and workflows
DISPLAY_INCLUDE_NESTED_SUBWORKFLOWS = True
DISPLAY_SUBWORKFLOWS_AS_PANEL = True
DISPLAY_INCLUDE_JOB_STATE_BREAKDOWN = True
DISPLAY_DIVIDER = "◆"
# bar.* are Rich defaults for these values.
DISPLAY_STYLE_BAR_BACK: StyleType = "bar.back"
DISPLAY_STYLE_BAR_FINISHED: StyleType = "bar.finished"
DISPLAY_STYLE_BAR_COMPLETE: StyleType = "bar.complete"

DISPLAY_STYLE_INITIALIZING: StyleType = "cyan"
DISPLAY_STYLE_OK: StyleType = "green"
DISPLAY_STYLE_RUNNING: StyleType = "green"
DISPLAY_STYLE_ERROR: StyleType = "red"

# Rich default style - a magenta
DISPLAY_STYLE_PERCENT: StyleType = "progress.percentage"

DISPLAY_STYLE_HEADER: StyleType = "bold"
DISPLAY_STYLE_SUBWORKFLOW_HEADER: StyleType = "bold"

DISPLAY_LABEL_HEADER_PREFIX = "Invocation "
DISPLAY_LABEL_SUBWORKFLOW_HEADER_PREFIX = "Subworkflow Invocation "
DISPLAY_LABEL_PROGRESS_STEPS = "Steps"
DISPLAY_LABEL_PROGRESS_JOBS = "Jobs"
DISPLAY_LABEL_PROGRESS_SUBWORKFLOWS = "SubWFs"
DISPLAY_LABEL_JOB_STATES_PREFIX = "Job States"

DISPLAY_ICON_STATE_OK = "🟢"
DISPLAY_ICON_STATE_ERRORS = "🔴"
DISPLAY_ICON_STATE_NEW = "🆕"
DISPLAY_ICON_STATE_QUEUED = "⏳"
DISPLAY_ICON_STATE_RUNNING = "👟"
DISPLAY_ICON_STATE_PAUSED = "⏸️"


class DisplayConfiguration(BaseModel):
    include_nested_subworkflows: bool = DISPLAY_INCLUDE_NESTED_SUBWORKFLOWS
    include_job_state_breakdown: bool = DISPLAY_INCLUDE_JOB_STATE_BREAKDOWN
    subworkflows_as_panel: bool = DISPLAY_SUBWORKFLOWS_AS_PANEL
    divider: str = DISPLAY_DIVIDER
    style_bar_back: StyleType = DISPLAY_STYLE_BAR_BACK
    style_bar_complete: StyleType = DISPLAY_STYLE_BAR_COMPLETE
    style_bar_finished: StyleType = DISPLAY_STYLE_BAR_FINISHED

    style_percent: StyleType = DISPLAY_STYLE_PERCENT

    style_initializing: StyleType = DISPLAY_STYLE_INITIALIZING
    style_ok: StyleType = DISPLAY_STYLE_OK
    style_running: StyleType = DISPLAY_STYLE_RUNNING
    style_error: StyleType = DISPLAY_STYLE_ERROR

    style_header: StyleType = DISPLAY_STYLE_HEADER
    style_subworkflow_header: StyleType = DISPLAY_STYLE_SUBWORKFLOW_HEADER

    label_header_prefix: str = DISPLAY_LABEL_HEADER_PREFIX
    label_subworkflow_header_prefix: str = DISPLAY_LABEL_SUBWORKFLOW_HEADER_PREFIX
    label_progress_steps: str = DISPLAY_LABEL_PROGRESS_STEPS
    label_progress_jobs: str = DISPLAY_LABEL_PROGRESS_JOBS
    label_progress_subworkflows: str = DISPLAY_LABEL_PROGRESS_SUBWORKFLOWS
    label_job_states_prefix: str = DISPLAY_LABEL_JOB_STATES_PREFIX

    icon_state_ok: str = DISPLAY_ICON_STATE_OK
    icon_state_errors: str = DISPLAY_ICON_STATE_ERRORS
    icon_state_new: str = DISPLAY_ICON_STATE_NEW
    icon_state_queued: str = DISPLAY_ICON_STATE_QUEUED
    icon_state_running: str = DISPLAY_ICON_STATE_RUNNING
    icon_state_paused: str = DISPLAY_ICON_STATE_PAUSED
