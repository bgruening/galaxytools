"""
Contains possible interactions with the Galaxy workflow invocations
"""

import logging
from typing import (
    Any,
    Optional,
    TYPE_CHECKING,
)

import requests

from bioblend import (
    CHUNK_SIZE,
    ConnectionError,
    NotReady,
    wait_on,
)
from bioblend.galaxy.client import Client
from bioblend.galaxy.workflows import InputsBy

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

log = logging.getLogger(__name__)

INVOCATION_TERMINAL_STATES = {"cancelled", "failed", "scheduled"}
# Invocation non-terminal states are: "cancelling", "new", "ready"


class InvocationClient(Client):
    gi: "GalaxyInstance"
    module = "invocations"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_invocations(
        self,
        workflow_id: Optional[str] = None,
        history_id: Optional[str] = None,
        user_id: Optional[str] = None,
        include_terminal: bool = True,
        limit: Optional[int] = None,
        *,
        offset: Optional[int] = None,
        view: str = "collection",
        step_details: bool = False,
        job_id: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """
        Get all workflow invocations, or select a subset by specifying optional
        arguments for filtering (e.g. a workflow ID).

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID to filter on

        :type history_id: str
        :param history_id: Encoded history ID to filter on

        :type user_id: str
        :param user_id: Encoded user ID to filter on. This must be
                        your own user ID if your are not an admin user.

        :type include_terminal: bool
        :param include_terminal: Whether to include invocations in terminal
          state.

        :type limit: int
        :param limit: Maximum number of invocations to return.

        :type offset: int
        :param offset: Number of invocations to skip. Return invocations
           starting from item offset+1.

        :type view: str
        :param view: Level of detail to return per invocation, either
                     'element' or 'collection'.

        :type step_details: bool
        :param step_details: If 'view' is 'element', also include details
                             on individual steps.

        :type job_id: str
        :param job_id: Encoded job ID to filter on.

        :rtype: list
        :return: A list of workflow invocations.
          For example::

            [{'history_id': '2f94e8ae9edff68a',
              'id': 'df7a1f0c02a5b08e',
              'model_class': 'WorkflowInvocation',
              'state': 'new',
              'update_time': '2015-10-31T22:00:22',
              'uuid': 'c8aa2b1c-801a-11e5-a9e5-8ca98228593c',
              'workflow_id': '03501d7626bd192f'}]
        """
        params = {"include_terminal": include_terminal, "view": view, "step_details": step_details}
        if workflow_id:
            params["workflow_id"] = workflow_id
        if history_id:
            params["history_id"] = history_id
        if job_id:
            params["job_id"] = job_id
        if user_id:
            params["user_id"] = user_id
        if limit is not None:
            params["limit"] = limit
        if offset is not None:
            params["offset"] = offset
        return self._get(params=params)

    def show_invocation(self, invocation_id: str) -> dict[str, Any]:
        """
        Get a workflow invocation dictionary representing the scheduling of a
        workflow. This dictionary may be sparse at first (missing inputs and
        invocation steps) and will become more populated as the workflow is
        actually scheduled.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: dict
        :return: The workflow invocation.
          For example::

            {'history_id': '2f94e8ae9edff68a',
             'id': 'df7a1f0c02a5b08e',
             'inputs': {'0': {'id': 'a7db2fac67043c7e',
               'src': 'hda',
               'uuid': '7932ffe0-2340-4952-8857-dbaa50f1f46a'}},
             'model_class': 'WorkflowInvocation',
             'state': 'ready',
             'steps': [{'action': None,
               'id': 'd413a19dec13d11e',
               'job_id': None,
               'model_class': 'WorkflowInvocationStep',
               'order_index': 0,
               'state': None,
               'update_time': '2015-10-31T22:00:26',
               'workflow_step_id': 'cbbbf59e8f08c98c',
               'workflow_step_label': None,
               'workflow_step_uuid': 'b81250fd-3278-4e6a-b269-56a1f01ef485'},
              {'action': None,
               'id': '2f94e8ae9edff68a',
               'job_id': 'e89067bb68bee7a0',
               'model_class': 'WorkflowInvocationStep',
               'order_index': 1,
               'state': 'new',
               'update_time': '2015-10-31T22:00:26',
               'workflow_step_id': '964b37715ec9bd22',
               'workflow_step_label': None,
               'workflow_step_uuid': 'e62440b8-e911-408b-b124-e05435d3125e'}],
             'update_time': '2015-10-31T22:00:26',
             'uuid': 'c8aa2b1c-801a-11e5-a9e5-8ca98228593c',
             'workflow_id': '03501d7626bd192f'}
        """
        url = self._make_url(invocation_id)
        return self._get(url=url)

    def rerun_invocation(
        self,
        invocation_id: str,
        inputs_update: Optional[dict] = None,
        params_update: Optional[dict] = None,
        history_id: Optional[str] = None,
        history_name: Optional[str] = None,
        import_inputs_to_history: bool = False,
        replacement_params: Optional[dict] = None,
        allow_tool_state_corrections: bool = False,
        inputs_by: Optional[InputsBy] = None,
        parameters_normalized: bool = False,
        resource_params: Optional[dict[str, Any]] = None,
        use_cached_job: bool = True,
    ) -> dict[str, Any]:
        """
        Rerun a workflow invocation. For more extensive documentation of all
        parameters, see the ``gi.workflows.invoke_workflow()`` method.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID to be rerun

        :type inputs_update: dict
        :param inputs_update: If different inputs should be used to the original
          invocation, this should contain a mapping of workflow inputs to the new
          datasets and dataset collections. Watch out for conflict with the
          legacy params_update.

        :type params_update: dict
        :param params_update: If different non-dataset tool parameters should be
          used to the original invocation, this should contain a mapping of the
          new parameter values.

        :type history_id: str
        :param history_id: The encoded history ID where to store the workflow
          outputs. Alternatively, ``history_name`` may be specified to create a
          new history.

        :type history_name: str
        :param history_name: Create a new history with the given name to store
          the workflow outputs. If both ``history_id`` and ``history_name`` are
          provided, ``history_name`` is ignored. If neither is specified, a new
          'Unnamed history' is created.

        :type import_inputs_to_history: bool
        :param import_inputs_to_history: If ``True``, used workflow inputs will
          be imported into the history. If ``False``, only workflow outputs will
          be visible in the given history.

        :type allow_tool_state_corrections: bool
        :param allow_tool_state_corrections: If True, allow Galaxy to fill in
          missing tool state when running workflows. This may be useful for
          workflows using tools that have changed over time or for workflows
          built outside of Galaxy with only a subset of inputs defined.

        :type replacement_params: dict
        :param replacement_params: pattern-based replacements for post-job
          actions

        :type inputs_by: str
        :param inputs_by: Determines how inputs are referenced. Can be
          "step_index|step_uuid" (default), "step_index", "step_id", "step_uuid", or "name".

        :type parameters_normalized: bool
        :param parameters_normalized: Whether Galaxy should normalize the input
          parameters to ensure everything is referenced by a numeric step ID.
          Default is ``False``, but when setting parameters for a subworkflow,
          ``True`` is required.

        :type resource_params: dict
        :param resource_params: A dictionary containing the resource parameters
          to be used for this workflow run.

        :type use_cached_job: bool
        :param use_cached_job: Whether to use cached jobs for the workflow
          invocation.

        :rtype: dict
        :return: A dict describing the new workflow invocation.

        .. note::
          This method works only on Galaxy 21.01 or later.
        """
        try:
            payload = self.get_invocation_request(invocation_id)
        except ConnectionError as e:
            if e.status_code != 404:
                raise
            # Galaxy release_24.1 or earlier
            invocation = self.show_invocation(invocation_id)
            workflow_step_id_to_index = {
                step["workflow_step_id"]: index for index, step in enumerate(invocation["steps"])
            }
            # Merge input_step_parameters (indexed by label) into inputs (indexed by step index)
            inputs = invocation["inputs"]
            for param_input_dict in invocation["input_step_parameters"].values():
                workflow_step_id = param_input_dict["workflow_step_id"]
                workflow_step_index = workflow_step_id_to_index[workflow_step_id]
                inputs[str(workflow_step_index)] = param_input_dict
            payload = {
                "inputs": inputs,
                "instance": True,
                "workflow_id": invocation["workflow_id"],
            }
        else:
            # Drop history_id from the payload as we will set history later
            payload.pop("history_id")
        workflow_id = payload["workflow_id"]
        if inputs_update:
            if payload.get("inputs") is None:
                payload["inputs"] = {}
            payload["inputs"].update(inputs_update)
        if params_update:
            if payload.get("parameters") is None:
                payload["parameters"] = {}
            payload["parameters"].update(params_update)
        if replacement_params:
            payload["replacement_params"] = replacement_params
        if history_id:
            payload["history"] = f"hist_id={history_id}"
        elif history_name:
            payload["history"] = history_name
        if not import_inputs_to_history:
            payload["no_add_to_history"] = True
        if allow_tool_state_corrections:
            payload["allow_tool_state_corrections"] = allow_tool_state_corrections
        if inputs_by is not None:
            payload["inputs_by"] = inputs_by
        if parameters_normalized:
            payload["parameters_normalized"] = parameters_normalized
        if resource_params:
            payload["resource_params"] = resource_params
        payload["use_cached_job"] = use_cached_job
        url = "/".join((self.gi.url, "workflows", workflow_id, "invocations"))
        return self.gi.make_post_request(url=url, payload=payload)

    def cancel_invocation(self, invocation_id: str) -> dict[str, Any]:
        """
        Cancel the scheduling of a workflow.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: dict
        :return: The workflow invocation being cancelled
        """
        url = self._make_url(invocation_id)
        return self._delete(url=url)

    def show_invocation_step(self, invocation_id: str, step_id: str) -> dict[str, Any]:
        """
        See the details of a particular workflow invocation step.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :type step_id: str
        :param step_id: Encoded workflow invocation step ID

        :rtype: dict
        :return: The workflow invocation step.
          For example::

            {'action': None,
             'id': '63cd3858d057a6d1',
             'job_id': None,
             'model_class': 'WorkflowInvocationStep',
             'order_index': 2,
             'state': None,
             'update_time': '2015-10-31T22:11:14',
             'workflow_step_id': '52e496b945151ee8',
             'workflow_step_label': None,
             'workflow_step_uuid': '4060554c-1dd5-4287-9040-8b4f281cf9dc'}
        """
        url = self._invocation_step_url(invocation_id, step_id)
        return self._get(url=url)

    def run_invocation_step_action(self, invocation_id: str, step_id: str, action: Any) -> dict[str, Any]:
        """Execute an action for an active workflow invocation step. The
        nature of this action and what is expected will vary based on the
        the type of workflow step (the only currently valid action is True/False
        for pause steps).

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :type step_id: str
        :param step_id: Encoded workflow invocation step ID

        :type action: object
        :param action: Action to use when updating state, semantics depends on
           step type.

        :rtype: dict
        :return: Representation of the workflow invocation step
        """
        url = self._invocation_step_url(invocation_id, step_id)
        payload = {"action": action}
        return self._put(payload=payload, url=url)

    def get_invocation_summary(self, invocation_id: str) -> dict[str, Any]:
        """
        Get a summary of an invocation, stating the number of jobs which
        succeed, which are paused and which have errored.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: dict
        :return: The invocation summary.
          For example::

            {'states': {'paused': 4, 'error': 2, 'ok': 2},
             'model': 'WorkflowInvocation',
             'id': 'a799d38679e985db',
             'populated_state': 'ok'}
        """
        url = self._make_url(invocation_id) + "/jobs_summary"
        return self._get(url=url)

    def get_invocation_step_jobs_summary(self, invocation_id: str) -> list[dict[str, Any]]:
        """
        Get a detailed summary of an invocation, listing all jobs with
        their job IDs and current states.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: list of dicts
        :return: The invocation step jobs summary.
          For example::

            [{'id': 'e85a3be143d5905b',
              'model': 'Job',
              'populated_state': 'ok',
              'states': {'ok': 1}},
             {'id': 'c9468fdb6dc5c5f1',
              'model': 'Job',
              'populated_state': 'ok',
              'states': {'running': 1}},
             {'id': '2a56795cad3c7db3',
              'model': 'Job',
              'populated_state': 'ok',
              'states': {'new': 1}}]
        """
        url = self._make_url(invocation_id) + "/step_jobs_summary"
        return self._get(url=url)

    def get_invocation_request(self, invocation_id: str) -> dict[str, Any]:
        """
        Get a request dict for an invocation.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: dict
        :return: The invocation request.

        .. note::
          This method works only on Galaxy 24.2 or later.
        """
        url = self._make_url(invocation_id) + "/request"
        return self._get(url=url)

    def get_invocation_report(self, invocation_id: str) -> dict[str, Any]:
        """
        Get a Markdown report for an invocation.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: dict
        :return: The invocation report.
          For example::

            {'markdown': '\\n# Workflow Execution Summary of Example workflow\\n\\n
             ## Workflow Inputs\\n\\n\\n## Workflow Outputs\\n\\n\\n
             ## Workflow\\n```galaxy\\n
             workflow_display(workflow_id=f2db41e1fa331b3e)\\n```\\n',
             'render_format': 'markdown',
             'workflows': {'f2db41e1fa331b3e': {'name': 'Example workflow'}}}
        """
        url = self._make_url(invocation_id) + "/report"
        return self._get(url=url)

    def get_invocation_report_pdf(self, invocation_id: str, file_path: str, chunk_size: int = CHUNK_SIZE) -> None:
        """
        Get a PDF report for an invocation.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :type file_path: str
        :param file_path: Path to save the report
        """
        url = self._make_url(invocation_id) + "/report.pdf"
        r = self.gi.make_get_request(url, stream=True)
        if r.status_code != 200:
            raise Exception(
                "Failed to get the PDF report, the necessary dependencies may not be installed on the Galaxy server."
            )
        with open(file_path, "wb") as outf:
            for chunk in r.iter_content(chunk_size):
                outf.write(chunk)

    # TODO: Move to a new ``bioblend.galaxy.short_term_storage`` module
    def _wait_for_short_term_storage(
        self, storage_request_id: str, maxwait: float = 60, interval: float = 3, json: bool = True, stream: bool = False
    ) -> Any:
        """
        Wait until a short term storage request is ready

        :type storage_request_id: str
        :param storage_request_id: Storage request ID to wait for.

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the storage request
          to become ready. After this time, a ``TimeoutException`` will be
          raised.

        :type interval: float
        :param interval: Time (in seconds) to wait between 2 consecutive checks.

        :return: The decoded response if ``json`` is set to ``True``, otherwise
          the response object
        """
        url = f"{self.gi.url}/short_term_storage/{storage_request_id}"
        is_ready_url = f"{url}/ready"

        def check_and_get_short_term_storage() -> Any:
            if self._get(url=is_ready_url):
                return self._get(url=url, json=json, stream=stream)
            raise NotReady(f"Storage request {storage_request_id} is not ready")

        return wait_on(check_and_get_short_term_storage, maxwait=maxwait, interval=interval)

    def get_invocation_archive(
        self,
        invocation_id: str,
        model_store_format: str = "tar.gz",
        include_files: bool = True,
        include_deleted: bool = False,
        include_hidden: bool = False,
        bco_merge_history_metadata: bool = False,
        bco_override_environment_variables: Optional[dict[str, Any]] = None,
        bco_override_empirical_error: Optional[dict[str, Any]] = None,
        bco_override_algorithmic_error: Optional[dict[str, Any]] = None,
        bco_override_xref: Optional[dict[str, Any]] = None,
        maxwait: float = 1200,
    ) -> requests.Response:
        """
        Get an invocation as a compressed archive.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :type model_store_format: str
        :param model_store_format: format of model store to export. Currently
          supported formats: tar.gz, tar, bag.zip, bag.tar, bag.tgz, rocrate.zip
          and bco.json .

        :type include_files: bool
        :param include_files: include materialized files in export when available

        :type include_deleted: bool
        :param include_deleted: include file contents for deleted datasets (if include_files is True).

        :type include_hidden: bool
        :param include_hidden: include file contents for hidden datasets (if include_files is True).

        :type bco_merge_history_metadata: bool
        :param bco_merge_history_metadata: when reading tags/annotations to generate BCO object include history metadata.

        :type bco_override_environment_variables: dict[str, Any]
        :param bco_override_environment_variables: override environment variables for 'execution_domain' when generating BioCompute object.

        :type bco_override_empirical_error: dict[str, Any]
        :param bco_override_empirical_error: override empirical error for 'error domain' when generating BioCompute object.

        :type bco_override_algorithmic_error: dict[str, Any]
        :param bco_override_algorithmic_error: override algorithmic error for 'error domain' when generating BioCompute object.

        :type bco_override_xref: dict[str, Any]
        :param bco_override_xref: Override xref for 'description domain' when generating BioCompute object.

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the invocation
          object to become ready. After this time, a ``TimeoutException`` will
          be raised.

        :rtype: requests.Response
        :return: request.Response
        """
        payload = {
            "model_store_format": model_store_format,
            "include_files": include_files,
            "include_deleted": include_deleted,
            "include_hidden": include_hidden,
            "bco_merge_history_metadata": bco_merge_history_metadata,
        }
        if bco_override_environment_variables is not None:
            payload["bco_override_environment_variables"] = bco_override_environment_variables
        if bco_override_empirical_error is not None:
            payload["bco_override_empirical_error"] = bco_override_empirical_error
        if bco_override_algorithmic_error is not None:
            payload["bco_override_algorithmic_error"] = bco_override_algorithmic_error
        if bco_override_xref is not None:
            payload["bco_override_xref"] = bco_override_xref

        url = self._make_url(invocation_id) + "/prepare_store_download"
        psd = self._post(url=url, payload=payload)
        storage_request_id = psd["storage_request_id"]
        return self._wait_for_short_term_storage(storage_request_id, maxwait=maxwait, json=False, stream=True)

    def get_invocation_biocompute_object(self, invocation_id: str, maxwait: float = 1200) -> dict[str, Any]:
        """
        Get a BioCompute object for an invocation.

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the BioCompute
          object to become ready. After this time, a ``TimeoutException`` will
          be raised.

        :rtype: dict
        :return: The BioCompute object
        """
        url = self._make_url(invocation_id) + "/prepare_store_download"
        payload = {"model_store_format": "bco.json"}
        try:
            psd = self._post(url=url, payload=payload)
        except ConnectionError as e:
            if e.status_code not in (400, 404):
                raise
            # Galaxy release_22.05 and earlier
            url = self._make_url(invocation_id) + "/biocompute"
            return self._get(url=url)
        else:
            storage_request_id = psd["storage_request_id"]
            return self._wait_for_short_term_storage(storage_request_id, maxwait=maxwait)

    def wait_for_invocation(
        self, invocation_id: str, maxwait: float = 12000, interval: float = 3, check: bool = True
    ) -> dict[str, Any]:
        """
        Wait until an invocation is in a terminal state.

        :type invocation_id: str
        :param invocation_id: Invocation ID to wait for.

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the invocation state
          to become terminal. After this time, a ``TimeoutException`` will be
          raised.

        :type interval: float
        :param interval: Time (in seconds) to wait between 2 consecutive checks.

        :type check: bool
        :param check: Whether to check if the invocation terminal state is
          'scheduled'.

        :rtype: dict
        :return: Details of the workflow invocation.
        """

        def check_and_get_invocation() -> dict[str, Any]:
            invocation = self.gi.invocations.show_invocation(invocation_id)
            state = invocation["state"]
            if state in INVOCATION_TERMINAL_STATES:
                if check and state != "scheduled":
                    raise Exception(f"Invocation {invocation_id} is in terminal state {state}")
                return invocation
            raise NotReady(f"Invocation {invocation_id} is in non-terminal state {state}")

        return wait_on(check_and_get_invocation, maxwait=maxwait, interval=interval)

    def _invocation_step_url(self, invocation_id: str, step_id: str) -> str:
        return "/".join((self._make_url(invocation_id), "steps", step_id))


__all__ = ("InvocationClient",)
