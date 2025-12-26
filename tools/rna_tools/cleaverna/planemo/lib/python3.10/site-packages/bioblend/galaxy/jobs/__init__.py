"""
Contains possible interactions with the Galaxy Jobs
"""

import logging
from typing import (
    Any,
    Literal,
    Optional,
    TYPE_CHECKING,
)

from bioblend import (
    NotReady,
    wait_on,
)
from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

log = logging.getLogger(__name__)

JOB_TERMINAL_STATES = {"deleted", "deleting", "error", "ok"}
# Job non-terminal states are: 'deleted_new', 'failed', 'new', 'paused',
# 'queued', 'resubmitted', 'running', 'upload', 'waiting'


class JobsClient(Client):
    module = "jobs"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_jobs(
        self,
        state: Optional[str] = None,
        history_id: Optional[str] = None,
        invocation_id: Optional[str] = None,
        tool_id: Optional[str] = None,
        workflow_id: Optional[str] = None,
        user_id: Optional[str] = None,
        date_range_min: Optional[str] = None,
        date_range_max: Optional[str] = None,
        limit: int = 500,
        offset: int = 0,
        user_details: bool = False,
        order_by: Literal["create_time", "update_time"] = "update_time",
    ) -> list[dict[str, Any]]:
        """
        Get all jobs, or select a subset by specifying optional arguments for
        filtering (e.g. a state).

        If the user is an admin, this will return jobs for all the users,
        otherwise only for the current user.

        :type state: str or list of str
        :param state: Job states to filter on.

        :type history_id: str
        :param history_id: Encoded history ID to filter on.

        :type invocation_id: string
        :param invocation_id: Encoded workflow invocation ID to filter on.

        :type tool_id: str or list of str
        :param tool_id: Tool IDs to filter on.

        :type workflow_id: string
        :param workflow_id: Encoded workflow ID to filter on.

        :type user_id: str
        :param user_id: Encoded user ID to filter on. Only admin users can
          access the jobs of other users.

        :type date_range_min: str
        :param date_range_min: Mininum job update date (in YYYY-MM-DD format) to
          filter on.

        :type date_range_max: str
        :param date_range_max: Maximum job update date (in YYYY-MM-DD format) to
          filter on.

        :type limit: int
        :param limit: Maximum number of jobs to return.

        :type offset: int
        :param offset: Number of jobs to skip. Return jobs starting from item
          offset+1.

        :type user_details: bool
        :param user_details: If ``True`` and the user is an admin, add the user
          email to each returned job dictionary.

        :type order_by: str
        :param order_by: Whether to order jobs by ``create_time`` or
          ``update_time`` (the default).

        :rtype: list of dict
        :return: Summary information for each selected job.
          For example::

            [{'create_time': '2014-03-01T16:16:48.640550',
              'exit_code': 0,
              'id': 'ebfb8f50c6abde6d',
              'model_class': 'Job',
              'state': 'ok',
              'tool_id': 'fasta2tab',
              'update_time': '2014-03-01T16:16:50.657399'},
             {'create_time': '2014-03-01T16:05:34.851246',
              'exit_code': 0,
              'id': '1cd8e2f6b131e891',
              'model_class': 'Job',
              'state': 'ok',
              'tool_id': 'upload1',
              'update_time': '2014-03-01T16:05:39.558458'}]

        .. note::
          The following parameters work only on Galaxy 21.05 or later: ``user_id``,
          ``limit``, ``offset``, ``workflow_id``, ``invocation_id``.
        """
        params: dict[str, Any] = {"limit": limit, "offset": offset}
        if state:
            params["state"] = state
        if history_id:
            params["history_id"] = history_id
        if invocation_id:
            params["invocation_id"] = invocation_id
        if tool_id:
            params["tool_id"] = tool_id
        if workflow_id:
            params["workflow_id"] = workflow_id
        if user_id:
            params["user_id"] = user_id
        if date_range_min:
            params["date_range_min"] = date_range_min
        if date_range_max:
            params["date_range_max"] = date_range_max
        if user_details:
            params["user_details"] = user_details
        if order_by:
            params["order_by"] = order_by
        return self._get(params=params)

    def show_job(self, job_id: str, full_details: bool = False) -> dict[str, Any]:
        """
        Get details of a given job of the current user.

        :type job_id: str
        :param job_id: job ID

        :type full_details: bool
        :param full_details: when ``True``, the complete list of details for the
          given job.

        :rtype: dict
        :return: A description of the given job.
          For example::

            {'create_time': '2014-03-01T16:17:29.828624',
             'exit_code': 0,
             'id': 'a799d38679e985db',
             'inputs': {'input': {'id': 'ebfb8f50c6abde6d', 'src': 'hda'}},
             'model_class': 'Job',
             'outputs': {'output': {'id': 'a799d38679e985db', 'src': 'hda'}},
             'params': {'chromInfo': '"/opt/galaxy-central/tool-data/shared/ucsc/chrom/?.len"',
                        'dbkey': '"?"',
                        'seq_col': '"2"',
                        'title_col': '["1"]'},
             'state': 'ok',
             'tool_id': 'tab2fasta',
             'update_time': '2014-03-01T16:17:31.930728'}
        """
        params = {}
        if full_details:
            params["full"] = full_details

        return self._get(id=job_id, params=params)

    def _build_for_rerun(self, job_id: str) -> dict[str, Any]:
        """
        Get details of a given job that can be used to rerun the corresponding tool.

        :type job_id: str
        :param job_id: job ID

        :rtype: dict
        :return: A description of the given job, with all parameters required to rerun.

        """
        url = "/".join((self._make_url(job_id), "build_for_rerun"))
        return self._get(url=url)

    def rerun_job(
        self,
        job_id: str,
        remap: bool = False,
        tool_inputs_update: Optional[dict[str, Any]] = None,
        history_id: Optional[str] = None,
    ) -> dict[str, Any]:
        """
        Rerun a job.

        :type job_id: str
        :param job_id: job ID

        :type remap: bool
        :param remap: when ``True``, the job output(s) will be remapped onto the dataset(s)
          created by the original job; if other jobs were waiting for this job to finish
          successfully, they will be resumed using the new outputs of this tool run. When
          ``False``, new job output(s) will be created. Note that if Galaxy does not permit
          remapping for the job in question, specifying ``True`` will result in an error.

        :type tool_inputs_update: dict
        :param tool_inputs_update: dictionary specifying any changes which should be
          made to tool parameters for the rerun job. This dictionary should have the same
          structure as is required when submitting the ``tool_inputs`` dictionary to
          ``gi.tools.run_tool()``, but only needs to include the inputs or parameters
          to be updated for the rerun job.

        :type history_id: str
        :param history_id: ID of the history in which the job should be executed; if
          not specified, the same history will be used as the original job run.

        :rtype: dict
        :return: Information about outputs and the rerun job

        .. note::
          This method works only on Galaxy 21.01 or later.
        """
        job_rerun_params = self._build_for_rerun(job_id)
        job_inputs = job_rerun_params["state_inputs"]

        if remap:
            if not job_rerun_params["job_remap"]:
                raise ValueError("remap was set to True, but this job is not remappable.")
            job_inputs["rerun_remap_job_id"] = job_id

        def update_inputs(inputs: dict[str, Any], tool_inputs_update: dict[str, Any]) -> None:
            """Recursively update inputs with tool_inputs_update"""
            for input_param, input_value in tool_inputs_update.items():
                if isinstance(input_value, dict):
                    update_inputs(inputs[input_param], input_value)
                else:
                    inputs[input_param] = input_value

        if tool_inputs_update:
            update_inputs(job_inputs, tool_inputs_update)

        url = "/".join((self.gi.url, "tools"))
        payload = {
            "history_id": history_id if history_id else job_rerun_params["history_id"],
            "tool_id": job_rerun_params["id"],
            "inputs": job_inputs,
            "input_format": "21.01",
        }
        return self._post(url=url, payload=payload)

    def get_state(self, job_id: str) -> str:
        """
        Display the current state for a given job of the current user.

        :type job_id: str
        :param job_id: job ID

        :rtype: str
        :return: state of the given job among the following values: `new`,
          `queued`, `running`, `waiting`, `ok`. If the state cannot be
          retrieved, an empty string is returned.

        .. versionadded:: 0.5.3
        """
        return self.show_job(job_id).get("state", "")

    def search_jobs(self, tool_id: str, inputs: dict[str, Any], state: Optional[str] = None) -> list[dict[str, Any]]:
        """
        Return jobs matching input parameters.

        :type tool_id: str
        :param tool_id: only return jobs associated with this tool ID

        :type inputs: dict
        :param inputs: return only jobs that have matching inputs

        :type state: str
        :param state: only return jobs in this state

        :rtype: list of dicts
        :return: Summary information for each matching job

        This method is designed to scan the list of previously run jobs and find
        records of jobs with identical input parameters and datasets. This can
        be used to minimize the amount of repeated work by simply recycling the
        old results.

        .. versionchanged:: 0.16.0
          Replaced the ``job_info`` parameter with separate ``tool_id``,
          ``inputs`` and ``state``.
        """
        job_info = {
            "tool_id": tool_id,
            "inputs": inputs,
        }
        if state:
            job_info["state"] = state
        url = self._make_url() + "/search"
        return self._post(url=url, payload=job_info)

    def get_metrics(self, job_id: str) -> list[dict[str, Any]]:
        """
        Return job metrics for a given job.

        :type job_id: str
        :param job_id: job ID

        :rtype: list
        :return: list containing job metrics

        .. note::
          Calling ``show_job()`` with ``full_details=True`` also returns the
          metrics for a job if the user is an admin. This method allows to fetch
          metrics even as a normal user as long as the Galaxy instance has the
          ``expose_potentially_sensitive_job_metrics`` option set to ``true`` in
          the ``config/galaxy.yml`` configuration file.
        """
        url = self._make_url(module_id=job_id) + "/metrics"
        return self._get(url=url)

    def cancel_job(self, job_id: str) -> bool:
        """
        Cancel a job, deleting output datasets.

        :type job_id: str
        :param job_id: job ID

        :rtype: bool
        :return: ``True`` if the job was successfully cancelled, ``False`` if
          it was already in a terminal state before the cancellation.
        """
        return self._delete(id=job_id)

    def report_error(self, job_id: str, dataset_id: str, message: str, email: Optional[str] = None) -> dict[str, Any]:
        """
        Report an error for a given job and dataset to the server administrators.

        :type job_id: str
        :param job_id: job ID

        :type dataset_id: str
        :param dataset_id: Dataset ID

        :type message: str
        :param message: Error message

        :type email: str
        :param email: Email for error report submission. If not specified, the email
          associated with the Galaxy user account is used by default.

        :rtype: dict
        :return: dict containing job error reply

        .. note::
          This method works only on Galaxy 20.01 or later.
        """
        payload = {
            "message": message,
            "dataset_id": dataset_id,
        }
        if email is not None:
            payload["email"] = email

        url = self._make_url(module_id=job_id) + "/error"
        return self._post(url=url, payload=payload)

    def get_common_problems(self, job_id: str) -> dict[str, Any]:
        """
        Query inputs and jobs for common potential problems that might
        have resulted in job failure.

        :type job_id: str
        :param job_id: job ID

        :rtype: dict
        :return: dict containing potential problems

        .. note::
          This method works only on Galaxy 19.05 or later.
        """
        url = self._make_url(module_id=job_id) + "/common_problems"
        return self._get(url=url)

    def get_inputs(self, job_id: str) -> list[dict[str, Any]]:
        """
        Get dataset inputs used by a job.

        :type job_id: str
        :param job_id: job ID

        :rtype: list of dicts
        :return: Inputs for the given job
        """
        url = self._make_url(module_id=job_id) + "/inputs"
        return self._get(url=url)

    def get_outputs(self, job_id: str) -> list[dict[str, Any]]:
        """
        Get dataset outputs produced by a job.

        :type job_id: str
        :param job_id: job ID

        :rtype: list of dicts
        :return: Outputs of the given job
        """
        url = self._make_url(module_id=job_id) + "/outputs"
        return self._get(url=url)

    def resume_job(self, job_id: str) -> list[dict[str, Any]]:
        """
        Resume a job if it is paused.

        :type job_id: str
        :param job_id: job ID

        :rtype: list of dicts
        :return: list of dictionaries containing output dataset associations
        """
        url = self._make_url(module_id=job_id) + "/resume"
        return self._put(url=url)

    def get_destination_params(self, job_id: str) -> dict[str, Any]:
        """
        Get destination parameters for a job, describing
        the environment and location where the job is run.

        :type job_id: str
        :param job_id: job ID

        :rtype: dict
        :return: Destination parameters for the given job

        .. note::
          This method works only on Galaxy 20.05 or later and if the user is a
          Galaxy admin.
        """
        url = self._make_url(module_id=job_id) + "/destination_params"
        return self._get(url=url)

    def show_job_lock(self) -> bool:
        """
        Show whether the job lock is active or not. If it is active,
        no jobs will dispatch on the Galaxy server.

        :rtype: bool
        :return: Status of the job lock

        .. note::
          This method works only on Galaxy 20.05 or later and if the user is a
          Galaxy admin.
        """
        url = self.gi.url + "/job_lock"
        response = self._get(url=url)
        return response["active"]

    def update_job_lock(self, active: bool = False) -> bool:
        """
        Update the job lock status by setting ``active`` to either
        ``True`` or ``False``. If ``True``, all job dispatching will
        be blocked.

        :rtype: bool
        :return: Updated status of the job lock

        .. note::
          This method works only on Galaxy 20.05 or later and if the user is a
          Galaxy admin.
        """
        payload = {
            "active": active,
        }
        url = self.gi.url + "/job_lock"
        response = self._put(url=url, payload=payload)
        return response["active"]

    def wait_for_job(
        self, job_id: str, maxwait: float = 12000, interval: float = 3, check: bool = True
    ) -> dict[str, Any]:
        """
        Wait until a job is in a terminal state.

        :type job_id: str
        :param job_id: job ID

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the job state to
          become terminal. After this time, a ``TimeoutException`` will be
          raised.

        :type interval: float
        :param interval: Time (in seconds) to wait between 2 consecutive checks.

        :type check: bool
        :param check: Whether to check if the job terminal state is 'ok'.

        :rtype: dict
        :return: Details of the given job.
        """

        def check_and_get_job() -> dict[str, Any]:
            job = self.show_job(job_id)
            state = job["state"]
            if state in JOB_TERMINAL_STATES:
                if check and state != "ok":
                    raise Exception(f"Job {job_id} is in terminal state {state}")
                return job
            raise NotReady(f"Job {job_id} is in non-terminal state {state}")

        return wait_on(check_and_get_job, maxwait=maxwait, interval=interval)
