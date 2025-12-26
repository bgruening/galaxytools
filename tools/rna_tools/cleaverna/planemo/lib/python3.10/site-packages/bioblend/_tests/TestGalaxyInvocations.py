import contextlib
import os
import tempfile
import time
import zipfile
from typing import Any

from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyInvocations(GalaxyTestBase.GalaxyTestBase):
    workflow_id: str
    pause_workflow_id: str
    x_random_lines_workflow_id: str

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        cls.workflow_id = cls.gi.workflows.import_workflow_from_local_path(path)["id"]
        path = test_util.get_abspath(os.path.join("data", "test_workflow_pause.ga"))
        cls.pause_workflow_id = cls.gi.workflows.import_workflow_from_local_path(path)["id"]
        path = test_util.get_abspath(os.path.join("data", "select_x_random_lines.ga"))
        cls.x_random_lines_workflow_id = cls.gi.workflows.import_workflow_from_local_path(path)["id"]

    def setUp(self):
        super().setUp()
        self.history_id = self.gi.histories.create_history(name="TestGalaxyInvocations")["id"]
        self.dataset_id = self._test_dataset(self.history_id)

    def tearDown(self):
        self.gi.histories.delete_history(self.history_id, purge=True)

    @test_util.skip_unless_galaxy("release_19.09")
    def test_cancel_invocation(self):
        invocation = self._invoke_pause_workflow()
        invocation_id = invocation["id"]
        invocations = self.gi.invocations.get_invocations()
        assert invocation_id in [inv["id"] for inv in invocations]

        invocation = self.gi.invocations.show_invocation(invocation_id)
        assert invocation["state"] in ["new", "ready"]

        invocation = self.gi.invocations.cancel_invocation(invocation_id)
        assert invocation["state"] in ["cancelled", "cancelling"]

    @test_util.skip_unless_galaxy("release_20.01")
    def test_get_invocations(self):
        invoc1 = self._invoke_workflow()

        # Run another workflow on the same history
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        workflow2_id = self.gi.workflows.import_workflow_from_local_path(path)["id"]
        dataset = {"src": "hda", "id": self.dataset_id}
        invoc2 = self.gi.workflows.invoke_workflow(
            workflow2_id, history_id=self.history_id, inputs={"Input 1": dataset, "Input 2": dataset}, inputs_by="name"
        )

        # Run the second workflow on another history
        hist2_id = self.gi.histories.create_history("hist2")["id"]
        invoc3 = self.gi.workflows.invoke_workflow(
            workflow2_id, history_id=hist2_id, inputs={"Input 1": dataset, "Input 2": dataset}, inputs_by="name"
        )

        for invoc in (invoc1, invoc2, invoc3):
            self.gi.invocations.wait_for_invocation(invoc["id"])

        # Test filtering by workflow ID
        invocs = self.gi.invocations.get_invocations(workflow_id=workflow2_id)
        assert len(invocs) == 2
        for invoc in invocs:
            assert invoc["workflow_id"] == workflow2_id

        # Test filtering by history ID
        for hist_id, expected_invoc_num in {self.history_id: 2, hist2_id: 1}.items():
            invocs = self.gi.invocations.get_invocations(history_id=hist_id)
            assert len(invocs) == expected_invoc_num
            for invoc in invocs:
                assert invoc["history_id"] == hist_id

        # Test limiting
        limit_invocs = self.gi.invocations.get_invocations(limit=2)
        assert len(limit_invocs) == 2

        self.gi.histories.delete_history(hist2_id, purge=True)

    @test_util.skip_unless_galaxy("release_21.05")
    def test_get_invocations_filtering_by_job_id(self):
        invocation = self._invoke_workflow()
        self.gi.invocations.wait_for_invocation(invocation["id"])
        hist_contents = self.gi.histories.show_history(self.history_id, contents=True)
        hist_last_dataset_id = hist_contents[-1]["id"]
        hist_last_job_id = self.gi.datasets.show_dataset(hist_last_dataset_id)["creating_job"]
        invocs = self.gi.invocations.get_invocations(job_id=hist_last_job_id)
        assert len(invocs) == 1
        assert invocs[0]["id"] == invocation["id"]

    @test_util.skip_unless_galaxy("release_19.09")
    def test_get_invocation_report(self):
        invocation = self._invoke_workflow()

        invocation_id = invocation["id"]
        report = self.gi.invocations.get_invocation_report(invocation_id)
        assert "paste_columns" in report["markdown"]
        with contextlib.suppress(Exception):
            # This can fail if dependencies as weasyprint are not installed on the Galaxy server
            ret = self.gi.invocations.get_invocation_report_pdf(invocation_id, "report.pdf")
            assert ret is None

    @test_util.skip_unless_galaxy("release_23.0")
    def test_get_invocation_archive(self):
        invocation = self._invoke_workflow()
        self.gi.invocations.wait_for_invocation(invocation["id"])
        with tempfile.TemporaryDirectory() as folder:
            file = f"{folder}/temp.rocrate.zip"
            response = self.gi.invocations.get_invocation_archive(
                invocation_id=invocation["id"],
                model_store_format="rocrate.zip",
            )
            with open(file, "bw") as archive:
                for chunk in response.iter_content(chunk_size=8192):
                    archive.write(chunk)
            # Verify file is not empty
            assert zipfile.is_zipfile(file)

    @test_util.skip_unless_galaxy("release_20.09")
    def test_get_invocation_biocompute_object(self):
        invocation = self._invoke_workflow()

        self.gi.invocations.wait_for_invocation(invocation["id"])
        biocompute_object = self.gi.invocations.get_invocation_biocompute_object(invocation["id"])
        assert len(biocompute_object["description_domain"]["pipeline_steps"]) == 1

    @test_util.skip_unless_galaxy("release_19.09")
    def test_get_invocation_jobs_summary(self):
        invocation = self._invoke_workflow()
        self.gi.invocations.wait_for_invocation(invocation["id"])
        jobs_summary = self.gi.invocations.get_invocation_summary(invocation["id"])
        assert jobs_summary["populated_state"] == "ok"
        step_jobs_summary = self.gi.invocations.get_invocation_step_jobs_summary(invocation["id"])
        assert len(step_jobs_summary) == 1
        assert step_jobs_summary[0]["populated_state"] == "ok"

    @test_util.skip_unless_galaxy("release_19.09")
    @test_util.skip_unless_tool("cat1")
    @test_util.skip_unless_tool("cat")
    def test_workflow_scheduling(self):
        invocation = self._invoke_pause_workflow()
        invocation_id = invocation["id"]

        def invocation_steps_by_order_index() -> dict[int, dict[str, Any]]:
            invocation = self.gi.invocations.show_invocation(invocation_id)
            return {s["order_index"]: s for s in invocation["steps"]}

        for _ in range(20):
            if 2 in invocation_steps_by_order_index():
                break
            time.sleep(0.5)
        steps = invocation_steps_by_order_index()
        pause_step = steps[2]
        assert self.gi.invocations.show_invocation_step(invocation_id, pause_step["id"])["action"] is None
        self.gi.invocations.run_invocation_step_action(invocation_id, pause_step["id"], action=True)
        assert self.gi.invocations.show_invocation_step(invocation_id, pause_step["id"])["action"]
        self.gi.invocations.wait_for_invocation(invocation["id"])

    @test_util.skip_unless_galaxy("release_21.01")
    def test_rerun_invocation(self):
        invocation = self._invoke_workflow()
        self.gi.invocations.wait_for_invocation(invocation["id"])
        rerun_invocation = self.gi.invocations.rerun_invocation(invocation["id"], import_inputs_to_history=True)
        self.gi.invocations.wait_for_invocation(rerun_invocation["id"])
        history = self.gi.histories.show_history(rerun_invocation["history_id"], contents=True)
        assert len(history) == 3

    @test_util.skip_unless_galaxy("release_21.01")
    def test_rerun_invocation_with_input_params(self):
        threeline_dataset_id = self._test_dataset(self.history_id, contents="A\nB\nC")
        invocation = self._invoke_x_random_lines_workflow(threeline_dataset_id)
        self.gi.invocations.wait_for_invocation(invocation["id"])
        rerun_invocation = self.gi.invocations.rerun_invocation(invocation["id"], history_id=self.history_id)
        self.gi.invocations.wait_for_invocation(rerun_invocation["id"])

    @test_util.skip_unless_galaxy("release_24.2")
    def test_rerun_invocation_with_input_params_changed(self):
        threeline_dataset_id = self._test_dataset(self.history_id, contents="A\nB\nC")
        invocation = self._invoke_x_random_lines_workflow(threeline_dataset_id)
        self.gi.invocations.wait_for_invocation(invocation["id"])
        inputs_update = {"how_many": 1}
        rerun_invocation = self.gi.invocations.rerun_invocation(
            invocation["id"], inputs_update=inputs_update, history_id=self.history_id
        )
        self.gi.invocations.wait_for_invocation(rerun_invocation["id"])
        rerun_request = self.gi.invocations.get_invocation_request(rerun_invocation["id"])
        assert rerun_request["inputs"]["how_many"] == 1

    def _invoke_workflow(self) -> dict[str, Any]:
        dataset = {"src": "hda", "id": self.dataset_id}

        return self.gi.workflows.invoke_workflow(
            self.workflow_id,
            inputs={"Input 1": dataset, "Input 2": dataset},
            history_id=self.history_id,
            inputs_by="name",
        )

    def _invoke_pause_workflow(self) -> dict[str, Any]:
        return self.gi.workflows.invoke_workflow(
            self.pause_workflow_id,
            inputs={"0": {"src": "hda", "id": self.dataset_id}},
            history_id=self.history_id,
        )

    def _invoke_x_random_lines_workflow(self, dataset_id: str) -> dict[str, Any]:
        return self.gi.workflows.invoke_workflow(
            self.x_random_lines_workflow_id,
            inputs={"from_what": {"src": "hda", "id": dataset_id}, "how_many": 2},
            history_id=self.history_id,
            inputs_by="name",
        )
