import os
from datetime import (
    datetime,
    timedelta,
)
from operator import itemgetter
from typing import Literal

from bioblend.galaxy.tools.inputs import (
    dataset,
    inputs,
)
from . import (
    GalaxyTestBase,
    test_util,
)


class TestGalaxyJobs(GalaxyTestBase.GalaxyTestBase):
    def setUp(self):
        super().setUp()
        self.history_id = self.gi.histories.create_history(name="TestGalaxyJobs")["id"]
        self.dataset_contents = "line 1\nline 2\rline 3\r\nline 4"
        self.dataset_id = self._test_dataset(self.history_id, contents=self.dataset_contents)

    def tearDown(self):
        self.gi.histories.delete_history(self.history_id, purge=True)

    @test_util.skip_unless_tool("cat1")
    def test_wait_for_job(self):
        tool_inputs = inputs().set("input1", dataset(self.dataset_id))
        tool_output = self.gi.tools.run_tool(history_id=self.history_id, tool_id="cat1", tool_inputs=tool_inputs)
        job_id = tool_output["jobs"][0]["id"]
        job = self.gi.jobs.wait_for_job(job_id)
        assert job["state"] == "ok"

    @test_util.skip_unless_tool("random_lines1")
    def test_get_jobs(self):
        self._run_tool()
        self._run_tool()

        jobs = self.gi.jobs.get_jobs(tool_id="random_lines1", history_id=self.history_id)
        assert len(jobs) == 2
        jobs = self.gi.jobs.get_jobs(history_id=self.history_id, state="failed")
        assert len(jobs) == 0
        yesterday = datetime.today() - timedelta(days=1)
        jobs = self.gi.jobs.get_jobs(date_range_max=yesterday.strftime("%Y-%m-%d"), history_id=self.history_id)
        assert len(jobs) == 0
        tomorrow = datetime.today() + timedelta(days=1)
        jobs = self.gi.jobs.get_jobs(date_range_min=tomorrow.strftime("%Y-%m-%d"))
        assert len(jobs) == 0
        jobs = self.gi.jobs.get_jobs(date_range_min=datetime.today().strftime("%Y-%m-%d"), history_id=self.history_id)
        assert len(jobs) == 3

    @test_util.skip_unless_galaxy("release_21.05")
    def test_get_jobs_with_filtering(self):
        path = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
        workflow_id = self.gi.workflows.import_workflow_from_local_path(path)["id"]
        dataset = {"src": "hda", "id": self.dataset_id}
        invocation1 = self.gi.workflows.invoke_workflow(
            workflow_id,
            inputs={"Input 1": dataset, "Input 2": dataset},
            history_id=self.history_id,
            inputs_by="name",
        )
        invocation2 = self.gi.workflows.invoke_workflow(
            workflow_id,
            inputs={"Input 1": dataset, "Input 2": dataset},
            history_id=self.history_id,
            inputs_by="name",
        )
        self.gi.invocations.wait_for_invocation(invocation1["id"])
        self.gi.invocations.wait_for_invocation(invocation2["id"])

        all_jobs = self.gi.jobs.get_jobs(history_id=self.history_id, order_by="create_time")
        assert len(all_jobs) == 3
        job1_id = all_jobs[1]["id"]
        jobs = self.gi.jobs.get_jobs(history_id=self.history_id, limit=1, offset=1, order_by="create_time")
        assert len(jobs) == 1
        assert jobs[0]["id"] == job1_id
        jobs = self.gi.jobs.get_jobs(invocation_id=invocation1["id"])
        assert len(jobs) == 1
        job_id_inv = jobs[0]["id"]
        jobs = self.gi.jobs.get_jobs(workflow_id=workflow_id)
        assert len(jobs) == 2
        assert job_id_inv in [job["id"] for job in jobs]

    @test_util.skip_unless_galaxy("release_21.01")
    @test_util.skip_unless_tool("random_lines1")
    def test_run_and_rerun_random_lines(self):
        original_output = self._run_tool(input_format="21.01")
        original_job_id = original_output["jobs"][0]["id"]

        rerun_output = self.gi.jobs.rerun_job(original_job_id)
        original_output_content = self.gi.datasets.download_dataset(original_output["outputs"][0]["id"])
        rerun_output_content = self.gi.datasets.download_dataset(rerun_output["outputs"][0]["id"])
        assert rerun_output_content == original_output_content

    @test_util.skip_unless_galaxy("release_21.01")
    @test_util.skip_unless_tool("Show beginning1")
    def test_rerun_and_remap(self):
        path = test_util.get_abspath(os.path.join("data", "select_first.ga"))
        wf = self.gi.workflows.import_workflow_from_local_path(path)
        wf_inputs = {
            "0": {"src": "hda", "id": self.dataset_id},
            "1": "-1",
        }
        invocation_id = self.gi.workflows.invoke_workflow(wf["id"], inputs=wf_inputs, history_id=self.history_id)["id"]
        invocation = self.gi.invocations.wait_for_invocation(invocation_id)
        job_steps = [step for step in invocation["steps"] if step["job_id"]]
        job_steps.sort(key=itemgetter("order_index"))
        try:
            self.gi.jobs.wait_for_job(job_steps[0]["job_id"])
        except Exception:
            pass  # indicates the job failed as expected
        else:
            raise Exception("The job should have failed")

        history_contents = self.gi.histories.show_history(self.history_id, contents=True)
        assert len(history_contents) == 3
        assert history_contents[1]["state"] == "error"
        assert history_contents[2]["state"] == "paused"

        # resume the paused step job
        resumed_outputs = self.gi.jobs.resume_job(job_steps[-1]["job_id"])
        assert resumed_outputs[0]["name"] == "out_file1"
        # the following does not pass stably - the job goes back to paused too quickly
        # history_contents_resumed = self.gi.histories.show_history(self.history_id, contents=True)
        # assert history_contents_resumed[2]["state"] != "paused"

        # now rerun and remap with correct input param
        failed_job_id = self.gi.datasets.show_dataset(history_contents[1]["id"])["creating_job"]
        tool_inputs_update = {"lineNum": "1"}
        rerun_job = self.gi.jobs.rerun_job(failed_job_id, remap=True, tool_inputs_update=tool_inputs_update)
        new_job_id = rerun_job["jobs"][0]["id"]

        # Wait for the last dataset in the history to be unpaused and complete
        last_dataset = self.gi.histories.show_history(self.history_id, contents=True)[-1]
        last_job_id = self.gi.datasets.show_dataset(last_dataset["id"])["creating_job"]
        self.gi.jobs.wait_for_job(new_job_id)
        self.gi.jobs.resume_job(last_job_id)  # last_job can get stuck on paused - resume it in case
        self.gi.jobs.wait_for_job(last_job_id)
        assert last_dataset["hid"] == 3
        assert last_dataset["id"] == history_contents[2]["id"]
        self._wait_and_verify_dataset(last_dataset["id"], b"line 1\tline 1\n")

    @test_util.skip_unless_tool("random_lines1")
    def test_get_common_problems(self):
        job_id = self._run_tool()["jobs"][0]["id"]
        response = self.gi.jobs.get_common_problems(job_id)
        assert response == {"has_duplicate_inputs": False, "has_empty_inputs": True}

    @test_util.skip_unless_tool("random_lines1")
    def test_get_inputs(self):
        job_id = self._run_tool()["jobs"][0]["id"]
        response = self.gi.jobs.get_inputs(job_id)
        assert response == [{"name": "input", "dataset": {"src": "hda", "id": self.dataset_id}}]

    @test_util.skip_unless_tool("random_lines1")
    def test_get_outputs(self):
        output = self._run_tool()
        job_id, output_id = output["jobs"][0]["id"], output["outputs"][0]["id"]
        response = self.gi.jobs.get_outputs(job_id)
        assert response == [{"name": "out_file1", "dataset": {"src": "hda", "id": output_id}}]

    @test_util.skip_unless_galaxy("release_20.05")
    @test_util.skip_unless_tool("random_lines1")
    def test_get_destination_params(self):
        job_id = self._run_tool()["jobs"][0]["id"]
        # In Galaxy 20.05 and 20.09 we need to wait for the job, otherwise
        # `get_destination_params()` receives a 500 error code. Fixed upstream
        # in https://github.com/galaxyproject/galaxy/commit/3e7f03cd1f229b8c9421ade02002728a33e131d8
        self.gi.jobs.wait_for_job(job_id)
        response = self.gi.jobs.get_destination_params(job_id)
        assert "Runner" in response
        assert "Runner Job ID" in response
        assert "Handler" in response

    @test_util.skip_unless_tool("random_lines1")
    def test_search_jobs(self):
        job_id = self._run_tool()["jobs"][0]["id"]
        inputs = {
            "num_lines": "1",
            "input": {"src": "hda", "id": self.dataset_id},
            "seed_source|seed_source_selector": "set_seed",
            "seed_source|seed": "asdf",
        }
        response = self.gi.jobs.search_jobs("random_lines1", inputs)
        assert job_id in [job["id"] for job in response]

    @test_util.skip_unless_galaxy("release_20.01")
    @test_util.skip_unless_tool("random_lines1")
    def test_report_error(self):
        output = self._run_tool()
        job_id, output_id = output["jobs"][0]["id"], output["outputs"][0]["id"]
        response = self.gi.jobs.report_error(job_id, output_id, "Test error")
        # expected response when the Galaxy server does not have mail configured
        assert response == {
            "messages": [
                [
                    "An error occurred sending the report by email: Mail is not configured for this Galaxy instance",
                    "danger",
                ]
            ]
        }

    @test_util.skip_unless_galaxy("release_20.05")
    def test_show_job_lock(self):
        status = self.gi.jobs.show_job_lock()
        assert not status

    @test_util.skip_unless_galaxy("release_20.05")
    def test_update_job_lock(self):
        status = self.gi.jobs.update_job_lock(active=True)
        assert status
        status = self.gi.jobs.update_job_lock(active=False)
        assert not status

    def test_cancel_job(self):
        job_id = self._run_tool()["jobs"][0]["id"]
        self.gi.jobs.cancel_job(job_id)
        job = self.gi.jobs.wait_for_job(job_id, check=False)
        assert job["state"] in ("deleted", "deleting")

    def _run_tool(self, input_format: Literal["21.01", "legacy"] = "legacy") -> dict:
        return super()._run_random_lines1(self.history_id, self.dataset_id, input_format=input_format)
