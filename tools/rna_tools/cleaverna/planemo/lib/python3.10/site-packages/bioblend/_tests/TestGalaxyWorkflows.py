import json
import os
import shutil
import tempfile
import time
from typing import (
    Any,
    Literal,
    Optional,
)

import pytest

from bioblend import ConnectionError
from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyWorkflows(GalaxyTestBase.GalaxyTestBase):
    @test_util.skip_unless_tool("cat1")
    @test_util.skip_unless_tool("cat")
    def test_workflow_scheduling(self):
        path = test_util.get_abspath(os.path.join("data", "test_workflow_pause.ga"))
        workflow = self.gi.workflows.import_workflow_from_local_path(path)
        workflow_id = workflow["id"]
        history_id = self.gi.histories.create_history(name="TestWorkflowState")["id"]

        invocations = self.gi.workflows.get_invocations(workflow_id)
        assert len(invocations) == 0

        # Try invalid invocation (no input)
        with pytest.raises(ConnectionError):
            self.gi.workflows.invoke_workflow(workflow["id"])

        dataset1_id = self._test_dataset(history_id)
        invocation = self.gi.workflows.invoke_workflow(
            workflow["id"],
            inputs={"0": {"src": "hda", "id": dataset1_id}},
        )
        assert invocation["state"] == "new"
        invocation_id = invocation["id"]
        invocations = self.gi.workflows.get_invocations(workflow_id)
        assert len(invocations) == 1
        assert invocations[0]["id"] == invocation_id

        def invocation_steps_by_order_index() -> dict[int, dict[str, Any]]:
            invocation = self.gi.workflows.show_invocation(workflow_id, invocation_id)
            return {s["order_index"]: s for s in invocation["steps"]}

        for _ in range(20):
            if 2 in invocation_steps_by_order_index():
                break
            time.sleep(0.5)

        invocation = self.gi.workflows.show_invocation(workflow_id, invocation_id)
        assert invocation["state"] == "ready"

        steps = invocation_steps_by_order_index()
        pause_step = steps[2]
        assert self.gi.workflows.show_invocation_step(workflow_id, invocation_id, pause_step["id"])["action"] is None
        self.gi.workflows.run_invocation_step_action(workflow_id, invocation_id, pause_step["id"], action=True)
        assert self.gi.workflows.show_invocation_step(workflow_id, invocation_id, pause_step["id"])["action"]
        for _ in range(20):
            invocation = self.gi.workflows.show_invocation(workflow_id, invocation_id)
            if invocation["state"] == "scheduled":
                break

            time.sleep(0.5)

        invocation = self.gi.workflows.show_invocation(workflow_id, invocation_id)
        assert invocation["state"] == "scheduled"

    def test_invoke_workflow_parameters_normalized(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns_subworkflow.ga"))
        workflow_id = self.gi.workflows.import_workflow_from_local_path(path)["id"]
        history_id = self.gi.histories.create_history(name="TestWorkflowInvokeParametersNormalized")["id"]
        dataset_id = self._test_dataset(history_id)
        with pytest.raises(ConnectionError):
            self.gi.workflows.invoke_workflow(
                workflow_id, inputs={"0": {"src": "hda", "id": dataset_id}}, params={"1": {"1|2": "comma"}}
            )
        self.gi.workflows.invoke_workflow(
            workflow_id,
            inputs={"0": {"src": "hda", "id": dataset_id}},
            params={"1": {"1|2": "comma"}},
            parameters_normalized=True,
        )

    @test_util.skip_unless_galaxy("release_19.09")
    @test_util.skip_unless_tool("cat1")
    @test_util.skip_unless_tool("cat")
    def test_cancelling_workflow_scheduling(self):
        path = test_util.get_abspath(os.path.join("data", "test_workflow_pause.ga"))
        workflow_id = self.gi.workflows.import_workflow_from_local_path(path)["id"]
        history_id = self.gi.histories.create_history(name="TestWorkflowState")["id"]
        dataset1_id = self._test_dataset(history_id)

        invocations = self.gi.workflows.get_invocations(workflow_id)
        assert len(invocations) == 0

        invocation = self.gi.workflows.invoke_workflow(
            workflow_id,
            inputs={"0": {"src": "hda", "id": dataset1_id}},
        )
        invocation_id = invocation["id"]
        invocations = self.gi.workflows.get_invocations(workflow_id)
        assert invocation_id in [inv["id"] for inv in invocations]

        invocation = self.gi.workflows.show_invocation(workflow_id, invocation_id)
        assert invocation["state"] in ["new", "ready"]

        invocation = self.gi.workflows.cancel_invocation(workflow_id, invocation_id)
        assert invocation["state"] in ["cancelled", "cancelling"]

    def test_import_export_workflow_from_local_path(self):
        with pytest.raises(TypeError):
            self.gi.workflows.import_workflow_from_local_path(None)  # type: ignore[arg-type]
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        imported_wf = self.gi.workflows.import_workflow_from_local_path(path)
        assert isinstance(imported_wf, dict)
        assert imported_wf["name"] == "paste_columns"
        assert imported_wf["url"].startswith("/api/workflows/")
        assert not imported_wf["deleted"]
        assert not imported_wf["published"]
        with pytest.raises(TypeError):
            self.gi.workflows.export_workflow_to_local_path(None, None, None)  # type: ignore[arg-type]
        export_dir = tempfile.mkdtemp(prefix="bioblend_test_")
        try:
            self.gi.workflows.export_workflow_to_local_path(imported_wf["id"], export_dir)
            dir_contents = os.listdir(export_dir)
            assert len(dir_contents) == 1
            export_path = os.path.join(export_dir, dir_contents[0])
            with open(export_path) as f:
                exported_wf_dict = json.load(f)
        finally:
            shutil.rmtree(export_dir)
        assert isinstance(exported_wf_dict, dict)

    def test_import_publish_workflow_from_local_path(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        imported_wf = self.gi.workflows.import_workflow_from_local_path(path, publish=True)
        assert isinstance(imported_wf, dict)
        assert not imported_wf["deleted"]
        assert imported_wf["published"]

    def test_import_other_users_published_workflow(self) -> None:
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        imported_wf = self.gi.workflows.import_workflow_from_local_path(path, publish=True)
        _, new_gi = test_util.new_user_gi(self.gi)
        imported_wf_by_new_user = new_gi.workflows.import_shared_workflow(imported_wf["id"])
        assert imported_wf_by_new_user["name"] == f"imported: {imported_wf['name']}"
        assert imported_wf_by_new_user["url"].startswith("/api/workflows/")
        assert not imported_wf_by_new_user["deleted"]
        assert not imported_wf_by_new_user["published"]

    def _import_export(self, style: Optional[Literal["ga", "format2"]] = None):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        with open(path) as f:
            wf_dict = json.load(f)
        imported_wf = self.gi.workflows.import_workflow_dict(wf_dict)
        assert isinstance(imported_wf, dict)
        assert imported_wf["name"] == "paste_columns"
        assert imported_wf["url"].startswith("/api/workflows/")
        assert not imported_wf["deleted"]
        assert not imported_wf["published"]
        exported_wf_dict = self.gi.workflows.export_workflow_dict(imported_wf["id"], style=style)
        assert isinstance(exported_wf_dict, dict)
        if style == "format2":
            assert exported_wf_dict["class"] == "GalaxyWorkflow"
        else:
            assert exported_wf_dict["a_galaxy_workflow"] == "true"

    def test_import_export_workflow_dict(self):
        self._import_export()

    def test_import_export_workflow_dict_ga(self):
        self._import_export("ga")

    def test_import_export_workflow_dict_format2(self):
        self._import_export("format2")

    def test_import_publish_workflow_dict(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        with open(path) as f:
            wf_dict = json.load(f)
        imported_wf = self.gi.workflows.import_workflow_dict(wf_dict, publish=True)
        assert isinstance(imported_wf, dict)
        assert not imported_wf["deleted"]
        assert imported_wf["published"]

    def test_get_workflows(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        workflow = self.gi.workflows.import_workflow_from_local_path(path)
        all_wfs = self.gi.workflows.get_workflows()
        assert len(all_wfs) > 0
        wfs_with_name = self.gi.workflows.get_workflows(name=workflow["name"])
        wf_list = [w for w in wfs_with_name if w["id"] == workflow["id"]]
        assert len(wf_list) == 1
        wf_data = wf_list[0]
        if "create_time" in workflow:  # Galaxy >= 20.01
            assert wf_data["create_time"] == workflow["create_time"]
        else:  # Galaxy <= 22.01
            assert wf_data["url"] == workflow["url"]

    def test_show_workflow(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        wf = self.gi.workflows.import_workflow_from_local_path(path)
        wf_data = self.gi.workflows.show_workflow(wf["id"])
        assert wf_data["id"] == wf["id"]
        assert wf_data["name"] == wf["name"]
        assert wf_data["url"] == wf["url"]
        assert wf_data["version"] == 0
        assert len(wf_data["steps"]) == 3
        assert wf_data["inputs"] is not None

    def test_update_workflow_name(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        wf = self.gi.workflows.import_workflow_from_local_path(path)
        new_name = "new name"
        updated_wf = self.gi.workflows.update_workflow(wf["id"], name=new_name)
        assert updated_wf["name"] == new_name

    @test_util.skip_unless_galaxy("release_21.01")
    def test_update_workflow_published(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        wf = self.gi.workflows.import_workflow_from_local_path(path)
        assert not wf["published"]
        updated_wf = self.gi.workflows.update_workflow(wf["id"], published=True)
        assert updated_wf["published"]
        updated_wf = self.gi.workflows.update_workflow(wf["id"], published=False)
        assert not updated_wf["published"]

    @test_util.skip_unless_galaxy("release_19.09")
    def test_extract_workflow_from_history(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        wf = self.gi.workflows.import_workflow_from_local_path(path)
        history_id = self.gi.histories.create_history(name="test_wf_invocation")["id"]
        dataset1_id = self._test_dataset(history_id)
        dataset = {"src": "hda", "id": dataset1_id}
        invocation_id = self.gi.workflows.invoke_workflow(
            wf["id"],
            inputs={"Input 1": dataset, "Input 2": dataset},
            history_id=history_id,
            inputs_by="name",
        )["id"]
        invocation = self.gi.invocations.wait_for_invocation(invocation_id)
        wf1 = self.gi.workflows.show_workflow(invocation["workflow_id"])
        datasets = self.gi.histories.show_history(invocation["history_id"], contents=True)
        dataset_hids = [dataset["hid"] for dataset in datasets]
        job_ids = [step["job_id"] for step in invocation["steps"] if step["job_id"]]

        for job_id in job_ids:
            self.gi.jobs.wait_for_job(job_id)

        new_workflow_name = "My new workflow!"
        wf2 = self.gi.workflows.extract_workflow_from_history(
            history_id=invocation["history_id"],
            workflow_name=new_workflow_name,
            job_ids=job_ids,
            dataset_hids=dataset_hids,
        )
        wf2 = self.gi.workflows.show_workflow(wf2["id"])
        assert wf2["name"] == new_workflow_name
        assert len(wf1["steps"]) == len(wf2["steps"])
        for i in range(len(wf1["steps"])):
            assert wf1["steps"][str(i)]["type"] == wf2["steps"][str(i)]["type"]
            assert wf1["steps"][str(i)]["tool_id"] == wf2["steps"][str(i)]["tool_id"]

    @test_util.skip_unless_galaxy("release_21.01")
    def test_refactor_workflow(self):
        actions: list[dict[str, Any]] = [
            {"action_type": "add_input", "type": "data", "label": "foo"},
            {"action_type": "update_step_label", "label": "bar", "step": {"label": "foo"}},
        ]
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        wf = self.gi.workflows.import_workflow_from_local_path(path)
        response = self.gi.workflows.refactor_workflow(wf["id"], actions, dry_run=True)
        assert len(response["action_executions"]) == len(actions)
        assert response["dry_run"] is True
        updated_steps = response["workflow"]["steps"]
        assert len(updated_steps) == 4
        assert {step["label"] for step in updated_steps.values()} == {"bar", None, "Input 1", "Input 2"}


class TestGalaxyWorkflowVersions(GalaxyTestBase.GalaxyTestBase):
    new_name: str
    workflow_id: str

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        cls.workflow_id = cls.gi.workflows.import_workflow_from_local_path(path)["id"]
        cls.new_name = "new name"
        cls.gi.workflows.update_workflow(cls.workflow_id, name=cls.new_name)

    @test_util.skip_unless_galaxy(
        "release_19.09"
    )  # due to Galaxy bug fixed in https://github.com/galaxyproject/galaxy/pull/9014
    def test_show_workflow_versions(self):
        updated_wf = self.gi.workflows.show_workflow(self.workflow_id)
        assert updated_wf["name"] == self.new_name
        assert updated_wf["version"] == 1
        wf_v0 = self.gi.workflows.show_workflow(self.workflow_id, version=0)
        assert wf_v0["name"] == "paste_columns"
        assert wf_v0["version"] == 0
        wf_v1 = self.gi.workflows.show_workflow(self.workflow_id, version=1)
        assert updated_wf == wf_v1

    def test_show_versions(self):
        versions = self.gi.workflows.show_versions(self.workflow_id)
        assert len(versions) == 2
        for i, version in enumerate(versions):
            assert version["version"] == i
            assert "update_time" in version
            assert "steps" in version

    @test_util.skip_unless_galaxy(
        "release_24.01"
    )  # due to Galaxy bug fixed in https://github.com/galaxyproject/galaxy/pull/18378
    def test_invoke_previous_version(self):
        history_id = self.gi.histories.create_history(name="test_wf_invocation")["id"]
        dataset1_id = self._test_dataset(history_id)
        dataset = {"src": "hda", "id": dataset1_id}
        invocation_id = self.gi.workflows.invoke_workflow(
            self.workflow_id,
            inputs={"Input 1": dataset, "Input 2": dataset},
            history_id=history_id,
            inputs_by="name",
            version=0,
        )["id"]
        self.gi.invocations.wait_for_invocation(invocation_id)
        # Try invalid invocation (wrong version)
        with pytest.raises(ConnectionError):
            self.gi.workflows.invoke_workflow(
                self.workflow_id,
                inputs={"Input 1": dataset, "Input 2": dataset},
                history_id=history_id,
                inputs_by="name",
                version=2,
            )
