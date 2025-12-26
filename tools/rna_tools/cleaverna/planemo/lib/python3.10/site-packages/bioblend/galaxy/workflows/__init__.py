"""
Contains possible interactions with the Galaxy Workflows
"""

import json
import os
from typing import (
    Any,
    Literal,
    Optional,
    TYPE_CHECKING,
)

import yaml

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

InputsBy = Literal["step_index|step_uuid", "step_index", "step_id", "step_uuid", "name"]


class WorkflowClient(Client):
    module = "workflows"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    # the 'deleted' option is not available for workflows
    def get_workflows(
        self, workflow_id: Optional[str] = None, name: Optional[str] = None, published: bool = False
    ) -> list[dict[str, Any]]:
        """
        Get all workflows, or select a subset by specifying optional arguments
        for filtering (e.g. a workflow name).

        :type name: str
        :param name: Workflow name to filter on.

        :type published: bool
        :param published: if ``True``, return also published workflows

        :rtype: list
        :return: A list of workflow dicts.
                 For example::

                   [{'id': '92c56938c2f9b315',
                     'name': 'Simple',
                     'url': '/api/workflows/92c56938c2f9b315'}]

        .. versionchanged:: 1.1.1
           Using the deprecated ``workflow_id`` parameter now raises a
           ``ValueError`` exception.
        """
        if workflow_id is not None:
            raise ValueError(
                "The workflow_id parameter has been removed, use the show_workflow() method to view details of a workflow for which you know the ID."
            )
        params: dict[str, Any] = {}
        if published:
            params["show_published"] = True
        workflows = self._get(params=params)
        if name is not None:
            workflows = [_ for _ in workflows if _["name"] == name]
        return workflows

    def show_workflow(
        self,
        workflow_id: str,
        version: Optional[int] = None,
        instance: Optional[bool] = None,
        legacy: Optional[bool] = None,
    ) -> dict[str, Any]:
        """
        Display information needed to run a workflow.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type version: int
        :param version: Workflow version to show

        :type instance: bool
        :param instance: treat ``workflow_id`` as a Workflow ID if True,
          otherwise treat it as a StoredWorkflow ID (the default). This
          parameter works only on Galaxy 20.01 or later.

        :type legacy: bool
        :param legacy: whether to use the legacy workflow format (default is
          False). Before Galaxy 24.0, passing False to this parameter was
          mistakenly equivalent to passing True.

        :rtype: dict
        :return: A description of the workflow and its inputs.
          For example::

            {'id': '92c56938c2f9b315',
             'inputs': {'23': {'label': 'Input Dataset', 'value': ''}},
             'name': 'Simple',
             'url': '/api/workflows/92c56938c2f9b315'}
        """
        params: dict[str, Any] = {}

        if version is not None:
            params["version"] = version
        if instance is not None:
            params["instance"] = instance
        if legacy is not None:
            params["legacy"] = legacy

        return self._get(id=workflow_id, params=params)

    def get_workflow_inputs(self, workflow_id: str, label: str) -> list[str]:
        """
        Get a list of workflow input IDs that match the given label.
        If no input matches the given label, an empty list is returned.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type label: str
        :param label: label to filter workflow inputs on

        :rtype: list
        :return: list of workflow inputs matching the label query
        """
        wf = self._get(id=workflow_id)
        inputs = wf["inputs"]
        return [id for id in inputs if inputs[id]["label"] == label]

    def import_workflow_dict(self, workflow_dict: dict[str, Any], publish: bool = False) -> dict[str, Any]:
        """
        Imports a new workflow given a dictionary representing a previously
        exported workflow.

        :type workflow_dict: dict
        :param workflow_dict: dictionary representing the workflow to be imported

        :type publish: bool
        :param publish:  if ``True`` the uploaded workflow will be published;
                         otherwise it will be visible only by the user which uploads it (default)

        :rtype: dict
        :return: Information about the imported workflow.
          For example::

            {'name': 'Training: 16S rRNA sequencing with mothur: main tutorial',
             'tags': [],
             'deleted': false,
             'latest_workflow_uuid': '368c6165-ccbe-4945-8a3c-d27982206d66',
             'url': '/api/workflows/94bac0a90086bdcf',
             'number_of_steps': 44,
             'published': false,
             'owner': 'jane-doe',
             'model_class': 'StoredWorkflow',
             'id': '94bac0a90086bdcf'}
        """
        payload = {"workflow": workflow_dict, "publish": publish}

        url = self._make_url() + "/upload"
        return self._post(url=url, payload=payload)

    def import_workflow_from_local_path(self, file_local_path: str, publish: bool = False) -> dict[str, Any]:
        """
        Imports a new workflow given the path to a file containing a previously
        exported workflow.

        :type file_local_path: str
        :param file_local_path: File to upload to the server for new workflow

        :type publish: bool
        :param publish:  if ``True`` the uploaded workflow will be published;
                         otherwise it will be visible only by the user which uploads it (default)

        :rtype: dict
        :return: Information about the imported workflow.
          For example::

            {'name': 'Training: 16S rRNA sequencing with mothur: main tutorial',
             'tags': [],
             'deleted': false,
             'latest_workflow_uuid': '368c6165-ccbe-4945-8a3c-d27982206d66',
             'url': '/api/workflows/94bac0a90086bdcf',
             'number_of_steps': 44,
             'published': false,
             'owner': 'jane-doe',
             'model_class': 'StoredWorkflow',
             'id': '94bac0a90086bdcf'}

        """
        with open(file_local_path) as fp:
            workflow_json = json.load(fp)

        return self.import_workflow_dict(workflow_json, publish)

    def import_shared_workflow(self, workflow_id: str) -> dict[str, Any]:
        """
        Imports a new workflow from the shared published workflows.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :rtype: dict
        :return: A description of the workflow.
          For example::

            {'id': 'ee0e2b4b696d9092',
             'model_class': 'StoredWorkflow',
             'name': 'Super workflow that solves everything!',
             'published': False,
             'tags': [],
             'url': '/api/workflows/ee0e2b4b696d9092'}
        """
        payload = {"shared_workflow_id": workflow_id}
        url = self._make_url()
        return self._post(url=url, payload=payload)

    def export_workflow_dict(
        self, workflow_id: str, version: Optional[int] = None, style: Optional[Literal["ga", "format2"]] = None
    ) -> dict[str, Any]:
        """
        Exports a workflow.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type version: int
        :param version: Workflow version to export

        :type style
        :param style: Either "ga" for the original JSON format or "format2" for
          the modern YAML gxformat2 format.

        :rtype: dict
        :return: Dictionary representing the requested workflow
        """
        params: dict[str, Any] = {}
        if version is not None:
            params["version"] = version
        if style:
            params["style"] = style
        url = "/".join((self._make_url(), "download", workflow_id))
        json = style != "format2"
        response = self._get(url=url, params=params, json=json)
        if not json:
            return yaml.safe_load(response.text)
        return response

    def export_workflow_to_local_path(
        self, workflow_id: str, file_local_path: str, use_default_filename: bool = True
    ) -> None:
        """
        Exports a workflow in JSON format to a given local path.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type file_local_path: str
        :param file_local_path: Local path to which the exported file will be saved.
                                (Should not contain filename if use_default_name=True)

        :type use_default_filename: bool
        :param use_default_filename: If the use_default_name parameter is True, the exported
          file will be saved as file_local_path/Galaxy-Workflow-%s.ga, where %s
          is the workflow name. If use_default_name is False, file_local_path
          is assumed to contain the full file path including filename.

        :rtype: None
        :return: None
        """
        workflow_dict = self.export_workflow_dict(workflow_id)

        if use_default_filename:
            filename = f"Galaxy-Workflow-{workflow_dict['name']}.ga"
            file_local_path = os.path.join(file_local_path, filename)

        with open(file_local_path, "w") as fp:
            json.dump(workflow_dict, fp)

    def update_workflow(self, workflow_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Update a given workflow.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type workflow: dict
        :param workflow: dictionary representing the workflow to be updated

        :type name: str
        :param name: New name of the workflow

        :type annotation: str
        :param annotation: New annotation for the workflow

        :type menu_entry: bool
        :param menu_entry: Whether the workflow should appear in the user's menu

        :type tags: list of str
        :param tags: Replace workflow tags with the given list

        :type published: bool
        :param published: Whether the workflow should be published or unpublished

        :rtype: dict
        :return: Dictionary representing the updated workflow
        """
        return self._put(payload=kwargs, id=workflow_id)

    def invoke_workflow(
        self,
        workflow_id: str,
        inputs: Optional[dict] = None,
        params: Optional[dict] = None,
        history_id: Optional[str] = None,
        history_name: Optional[str] = None,
        import_inputs_to_history: bool = False,
        replacement_params: Optional[dict] = None,
        allow_tool_state_corrections: bool = False,
        inputs_by: Optional[InputsBy] = None,
        parameters_normalized: bool = False,
        require_exact_tool_versions: bool = True,
        version: Optional[int] = None,
        use_cached_job: bool = False,
        parameters: Optional[dict] = None,
        instance: bool = False,
        resource_params: Optional[dict[str, Any]] = None,
        preferred_object_store_id: Optional[str] = None,
        preferred_intermediate_object_store_id: Optional[str] = None,
        preferred_outputs_object_store_id: Optional[str] = None,
    ) -> dict[str, Any]:
        """
        Invoke the workflow identified by ``workflow_id``. This will
        cause a workflow to be scheduled and return an object describing
        the workflow invocation.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type inputs: dict
        :param inputs: A mapping of workflow inputs to datasets and dataset collections.
                       The datasets source can be a LibraryDatasetDatasetAssociation (``ldda``),
                       LibraryDataset (``ld``), HistoryDatasetAssociation (``hda``), or
                       HistoryDatasetCollectionAssociation (``hdca``).

                       The map must be in the following format:
                       ``{'<input_index>': {'id': <encoded dataset ID>, 'src': '[ldda, ld, hda, hdca]'}}``
                       (e.g. ``{'2': {'id': '29beef4fadeed09f', 'src': 'hda'}}``)

                       This map may also be indexed by the UUIDs of the workflow steps,
                       as indicated by the ``uuid`` property of steps returned from the
                       Galaxy API. Alternatively workflow steps may be addressed by
                       the label that can be set in the workflow editor. If using
                       uuid or label you need to also set the ``inputs_by`` parameter
                       to ``step_uuid`` or ``name``.

        :type params: dict
        :param params: A mapping of non-datasets tool parameters (see below)
          A synonym to the "parameters" dict below. Both cannot be provided.

        :type history_id: str
        :param history_id: The encoded history ID where to store the workflow
          output. Alternatively, ``history_name`` may be specified to create a
          new history.

        :type history_name: str
        :param history_name: Create a new history with the given name to store
          the workflow output. If both ``history_id`` and ``history_name`` are
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
          actions (see below)

        :type inputs_by: str
        :param inputs_by: Determines how inputs are referenced. Can be
          "step_index|step_uuid" (default), "step_index", "step_id", "step_uuid", or "name".

        :type parameters_normalized: bool
        :param parameters_normalized: Whether Galaxy should normalize ``params``
          to ensure everything is referenced by a numeric step ID. Default is
          ``False``, but when setting ``params`` for a subworkflow, ``True`` is
          required.

        :type require_exact_tool_versions: bool
        :param require_exact_tool_versions: Whether invocation should fail if
          Galaxy does not have the exact tool versions. Default is ``True``.
          Parameter does not any effect for Galaxy versions < 22.05.

        :type version: int
        :param version: The version of the workflow to invoke. If omitted or
          None, the latest workflow version will be invoked.

        :type use_cached_job: bool
        :param use_cached_job: Whether to use cached jobs for the workflow
          invocation.

        :type parameters: dict
        :param parameters: A mapping of non-datasets tool parameters (see below)
          A synonym to the "params" dict above. Both cannot be provided.

        :type instance: bool
        :param instance: treat ``workflow_id`` as a Workflow ID if True,
          otherwise treat it as a StoredWorkflow ID (the default). This
          parameter works only on Galaxy 21.05 or later.

        :type resource_params: dict
        :param resource_params: A dictionary containing the resource parameters
          to be used for this workflow run.

        :type preferred_object_store_id: str
        :param preferred_object_store_id: The object store id where you want all
        outputs of this workflow run be stored.

        :type preferred_intermediate_object_store_id: str
        :param preferred_intermediate_object_store_id: The object store id where
        you want the intermediate outputs of this workflow run to be stored.
        Cannot be set if ``preferred_object_store_id`` is set.

        :type preferred_outputs_object_store_id: str
        :param preferred_outputs_object_store_id: The object store id where
        you want the priamry outputs of this workflow run to be stored.
        Cannot be set if ``preferred_object_store_id`` is set.

        :rtype: dict
        :return: A dict containing the workflow invocation describing the
          scheduling of the workflow. For example::

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

        The ``params`` dict should be specified as follows::

          {STEP_ID: PARAM_DICT, ...}

        where PARAM_DICT is::

          {PARAM_NAME: VALUE, ...}

        For backwards compatibility, the following (deprecated) format is
        also supported for ``params``::

          {TOOL_ID: PARAM_DICT, ...}

        in which case PARAM_DICT affects all steps with the given tool id.
        If both by-tool-id and by-step-id specifications are used, the
        latter takes precedence.

        Finally (again, for backwards compatibility), PARAM_DICT can also
        be specified as::

          {'param': PARAM_NAME, 'value': VALUE}

        Note that this format allows only one parameter to be set per step.

        For a ``repeat`` parameter, the names of the contained parameters needs
        to be specified as ``<repeat name>_<repeat index>|<param name>``, with
        the repeat index starting at 0. For example, if the tool XML contains::

          <repeat name="cutoff" title="Parameters used to filter cells" min="1">
              <param name="name" type="text" value="n_genes" label="Name of param...">
                  <option value="n_genes">n_genes</option>
                  <option value="n_counts">n_counts</option>
              </param>
              <param name="min" type="float" min="0" value="0" label="Min value"/>
          </repeat>

        then the PARAM_DICT should be something like::

          {...
           "cutoff_0|name": "n_genes",
           "cutoff_0|min": "2",
           "cutoff_1|name": "n_counts",
           "cutoff_1|min": "4",
           ...}

        At the time of this writing, it is not possible to change the number of
        times the contained parameters are repeated. Therefore, the parameter
        indexes can go from 0 to n-1, where n is the number of times the
        repeated element was added when the workflow was saved in the Galaxy UI.

        The ``replacement_params`` dict should map parameter names in
        post-job actions (PJAs) to their runtime values. For
        instance, if the final step has a PJA like the following::

          {'RenameDatasetActionout_file1': {'action_arguments': {'newname': '${output}'},
                                            'action_type': 'RenameDatasetAction',
                                            'output_name': 'out_file1'}}

        then the following renames the output dataset to 'foo'::

          replacement_params = {'output': 'foo'}

        see also `this email thread
        <http://lists.bx.psu.edu/pipermail/galaxy-dev/2011-September/006875.html>`_.

        .. warning::
          Historically, workflow invocation consumed a ``dataset_map``
          data structure that was indexed by unencoded workflow step IDs. These
          IDs would not be stable across Galaxy instances. The new ``inputs``
          property is instead indexed by either the ``order_index`` property
          (which is stable across workflow imports) or the step UUID which is
          also stable.
        """
        payload: dict[str, Any] = {
            "allow_tool_state_corrections": allow_tool_state_corrections,
            "require_exact_tool_versions": require_exact_tool_versions,
            "version": version,
            "use_cached_job": use_cached_job,
            "instance": instance,
        }

        if params and parameters:
            raise ValueError("You may specify either 'parameters' or 'params' but not both.")
        elif parameters:
            payload["parameters"] = parameters
        elif params:
            payload["parameters"] = params

        split_object_store_config = (
            preferred_outputs_object_store_id is not None or preferred_intermediate_object_store_id is not None
        )
        if split_object_store_config and preferred_object_store_id:
            raise ValueError(
                "You may specify either 'preferred_object_store_id' or one/both of 'preferred_outputs_object_store_id' and 'preferred_intermediate_object_store_id' but not both"
            )

        if preferred_object_store_id:
            payload["preferred_object_store_id"] = preferred_object_store_id
        if preferred_intermediate_object_store_id:
            payload["preferred_intermediate_object_store_id"] = preferred_intermediate_object_store_id
        if preferred_outputs_object_store_id:
            payload["preferred_outputs_object_store_id"] = preferred_outputs_object_store_id

        if inputs:
            payload["inputs"] = inputs
        if replacement_params:
            payload["replacement_params"] = replacement_params
        if history_id:
            payload["history"] = f"hist_id={history_id}"
        elif history_name:
            payload["history"] = history_name
        if not import_inputs_to_history:
            payload["no_add_to_history"] = True
        if inputs_by is not None:
            payload["inputs_by"] = inputs_by
        if parameters_normalized:
            payload["parameters_normalized"] = parameters_normalized
        if resource_params:
            payload["resource_params"] = resource_params

        url = self._invocations_url(workflow_id)
        return self._post(payload, url=url)

    def show_invocation(self, workflow_id: str, invocation_id: str) -> dict[str, Any]:
        """
        Get a workflow invocation object representing the scheduling of a
        workflow. This object may be sparse at first (missing inputs and
        invocation steps) and will become more populated as the workflow is
        actually scheduled.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

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
        url = self._invocation_url(workflow_id, invocation_id)
        return self._get(url=url)

    def get_invocations(self, workflow_id: str) -> list[dict[str, Any]]:
        """
        Get a list containing all the workflow invocations corresponding to the
        specified workflow.

        For more advanced filtering use InvocationClient.get_invocations().

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

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
        url = self._invocations_url(workflow_id)
        return self._get(url=url)

    def cancel_invocation(self, workflow_id: str, invocation_id: str) -> dict[str, Any]:
        """
        Cancel the scheduling of a workflow.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type invocation_id: str
        :param invocation_id: Encoded workflow invocation ID

        :rtype: dict
        :return: The workflow invocation being cancelled
        """
        url = self._invocation_url(workflow_id, invocation_id)
        return self._delete(url=url)

    def show_invocation_step(self, workflow_id: str, invocation_id: str, step_id: str) -> dict[str, Any]:
        """
        See the details of a particular workflow invocation step.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

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
        url = self._invocation_step_url(workflow_id, invocation_id, step_id)
        return self._get(url=url)

    def run_invocation_step_action(
        self, workflow_id: str, invocation_id: str, step_id: str, action: Any
    ) -> dict[str, Any]:
        """Execute an action for an active workflow invocation step. The
        nature of this action and what is expected will vary based on the
        the type of workflow step (the only currently valid action is True/False
        for pause steps).

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

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
        url = self._invocation_step_url(workflow_id, invocation_id, step_id)
        payload = {"action": action}
        return self._put(payload=payload, url=url)

    def delete_workflow(self, workflow_id: str) -> None:
        """
        Delete a workflow identified by `workflow_id`.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        .. warning::
            Deleting a workflow is irreversible in Galaxy versions < 23.01 - all
            workflow data will be permanently deleted.
        """
        self._delete(id=workflow_id)

    def refactor_workflow(
        self, workflow_id: str, actions: list[dict[str, Any]], dry_run: bool = False
    ) -> dict[str, Any]:
        """
        Refactor workflow with given actions.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :type actions: list of dicts
        :param actions: Actions to use for refactoring the workflow. The following
                        actions are supported: update_step_label, update_step_position,
                        update_output_label, update_name, update_annotation,
                        update_license, update_creator, update_report, add_step,
                        add_input, disconnect, connect, fill_defaults, fill_step_defaults,
                        extract_input, extract_legacy_parameter,
                        remove_unlabeled_workflow_outputs, upgrade_all_steps,
                        upgrade_subworkflow, upgrade_tool.

          An example value for the ``actions`` argument might be::

            actions = [
                {"action_type": "add_input", "type": "data", "label": "foo"},
                {"action_type": "update_step_label", "label": "bar", "step": {"label": "foo"}},
            ]

        :type dry_run: bool
        :param dry_run: When true, perform a dry run where the existing
                        workflow is preserved. The refactored workflow
                        is returned in the output of the method, but not saved
                        on the Galaxy server.

        :rtype: dict
        :return: Dictionary containing logged messages for the executed actions
                 and the refactored workflow.
        """
        payload = {
            "actions": actions,
            "dry_run": dry_run,
        }
        url = "/".join((self._make_url(workflow_id), "refactor"))
        return self._put(payload=payload, url=url)

    def extract_workflow_from_history(
        self,
        history_id: str,
        workflow_name: str,
        job_ids: Optional[list[str]] = None,
        dataset_hids: Optional[list[str]] = None,
        dataset_collection_hids: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Extract a workflow from a history.

        :type history_id: str
        :param history_id: Encoded history ID

        :type   workflow_name: str
        :param  workflow_name: Name of the workflow to create

        :type   job_ids: list
        :param  job_ids: Optional list of job IDs to filter the jobs to extract from the history

        :type   dataset_hids: list
        :param  dataset_hids: Optional list of dataset hids corresponding to workflow inputs
                             when extracting a workflow from history

        :type   dataset_collection_hids: list
        :param  dataset_collection_hids: Optional list of dataset collection hids corresponding to workflow inputs
                                        when extracting a workflow from history

        :rtype: dict
        :return: A description of the created workflow
        """
        payload = {
            "from_history_id": history_id,
            "job_ids": job_ids if job_ids else [],
            "dataset_ids": dataset_hids if dataset_hids else [],
            "dataset_collection_ids": dataset_collection_hids if dataset_collection_hids else [],
            "workflow_name": workflow_name,
        }
        return self._post(payload=payload)

    def show_versions(self, workflow_id: str) -> list[dict[str, Any]]:
        """
        Get versions for a workflow.

        :type workflow_id: str
        :param workflow_id: Encoded workflow ID

        :rtype: list of dicts
        :return: Ordered list of version descriptions for this workflow
        """
        url = self._make_url(workflow_id) + "/versions"
        return self._get(url=url)

    def _invocation_step_url(self, workflow_id: str, invocation_id: str, step_id: str) -> str:
        return "/".join((self._invocation_url(workflow_id, invocation_id), "steps", step_id))

    def _invocation_url(self, workflow_id: str, invocation_id: str) -> str:
        return "/".join((self._invocations_url(workflow_id), invocation_id))

    def _invocations_url(self, workflow_id: str) -> str:
        return "/".join((self._make_url(workflow_id), "invocations"))


__all__ = ("WorkflowClient",)
