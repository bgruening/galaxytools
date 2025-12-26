# pylint: disable=C0103,E1101
import json
import os
import shutil
import socket
import sys
import tarfile
import tempfile
import unittest
import uuid
from collections.abc import (
    Collection,
    Iterable,
)
from ssl import SSLError
from typing import (
    Any,
    Callable,
    Literal,
    Union,
)
from urllib.error import URLError
from urllib.request import urlopen

import pytest

import bioblend
from bioblend.galaxy import dataset_collections
from bioblend.galaxy.objects import (
    galaxy_instance,
    wrappers,
)
from . import test_util

bioblend.set_stream_logger("test", level="INFO")
socket.setdefaulttimeout(10.0)
SAMPLE_FN = test_util.get_abspath(os.path.join("data", "paste_columns.ga"))
SAMPLE_WF_COLL_FN = test_util.get_abspath(os.path.join("data", "paste_columns_collections.ga"))
SAMPLE_WF_PARAMETER_INPUT_FN = test_util.get_abspath(os.path.join("data", "workflow_with_parameter_input.ga"))
FOO_DATA = "foo\nbar\n"
FOO_DATA_2 = "foo2\nbar2\n"
SAMPLE_WF_DICT = {
    "deleted": False,
    "id": "9005c5112febe774",
    "inputs": {
        "571": {"label": "Input Dataset", "value": ""},
        "572": {"label": "Input Dataset", "value": ""},
    },
    "model_class": "StoredWorkflow",
    "name": "paste_columns",
    "owner": "user_foo",
    "published": False,
    "steps": {
        "571": {
            "id": 571,
            "input_steps": {},
            "tool_id": None,
            "tool_inputs": {"name": "Input Dataset"},
            "tool_version": None,
            "type": "data_input",
        },
        "572": {
            "id": 572,
            "input_steps": {},
            "tool_id": None,
            "tool_inputs": {"name": "Input Dataset"},
            "tool_version": None,
            "type": "data_input",
        },
        "573": {
            "id": 573,
            "input_steps": {
                "input1": {"source_step": 571, "step_output": "output"},
                "input2": {"source_step": 572, "step_output": "output"},
            },
            "tool_id": "Paste1",
            "tool_inputs": {
                "delimiter": '"T"',
                "input1": "null",
                "input2": "null",
            },
            "tool_version": "1.0.0",
            "type": "tool",
        },
    },
    "tags": [],
    "url": "/api/workflows/9005c5112febe774",
}
SAMPLE_INV_DICT: dict[str, Any] = {
    "history_id": "2f94e8ae9edff68a",
    "id": "df7a1f0c02a5b08e",
    "inputs": {"0": {"id": "a7db2fac67043c7e", "src": "hda", "uuid": "7932ffe0-2340-4952-8857-dbaa50f1f46a"}},
    "model_class": "WorkflowInvocation",
    "state": "ready",
    "steps": [
        {
            "action": None,
            "id": "d413a19dec13d11e",
            "job_id": None,
            "model_class": "WorkflowInvocationStep",
            "order_index": 0,
            "state": None,
            "update_time": "2015-10-31T22:00:26",
            "workflow_step_id": "cbbbf59e8f08c98c",
            "workflow_step_label": None,
            "workflow_step_uuid": "b81250fd-3278-4e6a-b269-56a1f01ef485",
        },
        {
            "action": None,
            "id": "2f94e8ae9edff68a",
            "job_id": "e89067bb68bee7a0",
            "model_class": "WorkflowInvocationStep",
            "order_index": 1,
            "state": "new",
            "update_time": "2015-10-31T22:00:26",
            "workflow_step_id": "964b37715ec9bd22",
            "workflow_step_label": None,
            "workflow_step_uuid": "e62440b8-e911-408b-b124-e05435d3125e",
        },
    ],
    "update_time": "2015-10-31T22:00:26",
    "uuid": "c8aa2b1c-801a-11e5-a9e5-8ca98228593c",
    "workflow_id": "03501d7626bd192f",
}


def is_reachable(url: str) -> bool:
    res = None
    try:
        res = urlopen(url, timeout=5)
    except (SSLError, URLError, socket.timeout):
        return False
    if res is not None:
        res.close()
    return True


def upload_from_fs(
    lib: wrappers.Library, bnames: Iterable[str], **kwargs: Any
) -> tuple[list[wrappers.LibraryDataset], list[str]]:
    tempdir = tempfile.mkdtemp(prefix="bioblend_test_")
    try:
        fnames = [os.path.join(tempdir, _) for _ in bnames]
        for fn in fnames:
            with open(fn, "w") as f:
                f.write(FOO_DATA)
        dss = lib.upload_from_galaxy_fs(fnames, **kwargs)
    finally:
        shutil.rmtree(tempdir)
    return dss, fnames


class MockWrapper(wrappers.Wrapper):
    BASE_ATTRS = ("a", "b")
    a: int
    b: list[int]

    def __init__(self, *args: Any, **kwargs: Any) -> None:
        super().__init__(*args, **kwargs)


class TestWrapper(unittest.TestCase):
    def setUp(self):
        self.d: dict[str, Any] = {"a": 1, "b": [2, 3], "c": {"x": 4}}
        with pytest.raises(TypeError):
            wrappers.Wrapper(self.d)
        self.w = MockWrapper(self.d)

    def test_initialize(self):
        for k in MockWrapper.BASE_ATTRS:
            assert getattr(self.w, k) == self.d[k]
        self.w.a = 222
        self.w.b[0] = 222
        assert self.w.a == 222
        assert self.w.b[0] == 222
        assert self.d["a"] == 1
        assert self.d["b"][0] == 2
        with pytest.raises(AttributeError):
            _ = self.w.foo  # type: ignore[attr-defined]
        with pytest.raises(AttributeError):
            self.w.foo = 0  # type: ignore[assignment]

    def test_taint(self):
        assert not self.w.is_modified
        self.w.a = 111  # pylint: disable=W0201
        assert self.w.is_modified

    def test_serialize(self):
        w = MockWrapper.from_json(self.w.to_json())
        assert w.wrapped == self.w.wrapped

    def test_clone(self):
        w = self.w.clone()
        assert w.wrapped == self.w.wrapped
        w.b[0] = 111
        assert self.w.b[0] == 2

    def test_kwargs(self):
        parent = MockWrapper({"a": 10})
        w = MockWrapper(self.d, parent=parent)
        assert w.parent is parent
        with pytest.raises(AttributeError):
            w.parent = 0  # type: ignore[assignment,misc]


@test_util.skip_unless_galaxy()
class GalaxyObjectsTestBase(unittest.TestCase):
    gi: galaxy_instance.GalaxyInstance

    @classmethod
    def setUpClass(cls) -> None:
        galaxy_key = os.environ["BIOBLEND_GALAXY_API_KEY"]
        galaxy_url = os.environ["BIOBLEND_GALAXY_URL"]
        cls.gi = galaxy_instance.GalaxyInstance(galaxy_url, api_key=galaxy_key)


class TestWorkflow(GalaxyObjectsTestBase):
    def setUp(self):
        self.wf = wrappers.Workflow(SAMPLE_WF_DICT)

    def test_initialize(self):
        assert self.wf.id == "9005c5112febe774"
        assert self.wf.name == "paste_columns"
        assert not self.wf.deleted
        assert self.wf.owner == "user_foo"
        assert not self.wf.published
        assert self.wf.tags == []
        assert self.wf.input_labels_to_ids == {"Input Dataset": {"571", "572"}}
        assert self.wf.tool_labels_to_ids == {"Paste1": {"573"}}
        assert self.wf.data_input_ids == {"571", "572"}
        assert self.wf.source_ids == {"571", "572"}
        assert self.wf.sink_ids == {"573"}

    def test_dag(self):
        inv_dag: dict[str, set[str]] = {}
        for h, tails in self.wf.dag.items():
            for t in tails:
                inv_dag.setdefault(str(t), set()).add(h)
        assert self.wf.inv_dag == inv_dag
        heads = set(self.wf.dag)
        assert heads == set.union(*self.wf.inv_dag.values())
        tails = set(self.wf.inv_dag)
        assert tails == set.union(*self.wf.dag.values())
        ids = self.wf.sorted_step_ids()
        assert set(ids) == heads | tails
        for h, tails in self.wf.dag.items():
            for t in tails:
                assert ids.index(h) < ids.index(t)

    def test_steps(self):
        steps = SAMPLE_WF_DICT["steps"]
        assert isinstance(steps, dict)
        for sid, s in self.wf.steps.items():
            assert isinstance(s, wrappers.Step)
            assert s.id == sid
            assert sid in steps
            assert s.parent is self.wf
        assert self.wf.data_input_ids == {"571", "572"}
        assert self.wf.tool_ids == {"573"}

    def test_taint(self):
        assert not self.wf.is_modified
        self.wf.steps["571"].tool_id = "foo"
        assert self.wf.is_modified

    def test_input_map(self):
        history = wrappers.History({}, gi=self.gi)
        library = wrappers.Library({}, gi=self.gi)
        hda = wrappers.HistoryDatasetAssociation({"id": "hda_id"}, container=history, gi=self.gi)
        ldda = wrappers.LibraryDatasetDatasetAssociation({"id": "ldda_id"}, container=library, gi=self.gi)
        input_map = self.wf._convert_input_map({"0": hda, "1": ldda, "2": {"id": "hda2_id", "src": "hda"}})
        assert input_map == {
            "0": {"id": "hda_id", "src": "hda"},
            "1": {"id": "ldda_id", "src": "ldda"},
            "2": {"id": "hda2_id", "src": "hda"},
        }


@test_util.skip_unless_galaxy("release_19.09")
class TestInvocation(GalaxyObjectsTestBase):
    dataset: wrappers.HistoryDatasetAssociation
    history: wrappers.History
    inv: wrappers.Invocation
    workflow: wrappers.Workflow
    workflow_pause: wrappers.Workflow

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        cls.inv = wrappers.Invocation(SAMPLE_INV_DICT, gi=cls.gi)
        with open(SAMPLE_FN) as f:
            cls.workflow = cls.gi.workflows.import_new(f.read())
        path_pause = test_util.get_abspath(os.path.join("data", "test_workflow_pause.ga"))
        with open(path_pause) as f:
            cls.workflow_pause = cls.gi.workflows.import_new(f.read())
        cls.history = cls.gi.histories.create(name="TestInvocation")
        cls.dataset = cls.history.paste_content("1\t2\t3")

    @classmethod
    def tearDownClass(cls):
        cls.history.delete(purge=True)

    def test_initialize(self):
        assert self.inv.workflow_id == "03501d7626bd192f"
        assert self.inv.history_id == "2f94e8ae9edff68a"
        assert self.inv.id == "df7a1f0c02a5b08e"
        assert self.inv.state == "ready"
        assert self.inv.update_time == "2015-10-31T22:00:26"
        assert self.inv.uuid == "c8aa2b1c-801a-11e5-a9e5-8ca98228593c"

    def test_initialize_steps(self):
        for step, step_dict in zip(self.inv.steps, SAMPLE_INV_DICT["steps"]):
            assert isinstance(step_dict, dict)
            assert isinstance(step, wrappers.InvocationStep)
            assert step.parent is self.inv
            assert step.id == step_dict["id"]
            assert step.job_id == step_dict["job_id"]
            assert step.order_index == step_dict["order_index"]
            assert step.state == step_dict["state"]
            assert step.update_time == step_dict["update_time"]
            assert step.workflow_step_id == step_dict["workflow_step_id"]
            assert step.workflow_step_label == step_dict["workflow_step_label"]
            assert step.workflow_step_uuid == step_dict["workflow_step_uuid"]

    def test_initialize_inputs(self):
        for i, input in enumerate(self.inv.inputs):
            assert input == {**SAMPLE_INV_DICT["inputs"][str(i)], "label": str(i)}

    def test_sorted_step_ids(self):
        assert self.inv.sorted_step_ids() == ["d413a19dec13d11e", "2f94e8ae9edff68a"]

    def test_step_states(self):
        assert self.inv.step_states() == {None, "new"}

    def test_number_of_steps(self):
        assert self.inv.number_of_steps() == 2

    def test_sorted_steps_by(self):
        assert len(self.inv.sorted_steps_by()) == 2
        steps = self.inv.sorted_steps_by(step_ids={"2f94e8ae9edff68a"})
        assert len(steps) == 1
        assert steps[0].id == "2f94e8ae9edff68a"
        assert self.inv.sorted_steps_by(step_ids={"unmatched_id"}) == []
        steps = self.inv.sorted_steps_by(states={"new"})
        assert len(steps) == 1
        assert steps[0].state == "new"
        assert self.inv.sorted_steps_by(states={"unmatched_state"}) == []
        steps = self.inv.sorted_steps_by(indices={0}, states={None, "new"})
        assert len(steps) == 1
        assert steps[0].order_index == 0
        assert self.inv.sorted_steps_by(indices={2}) == []

    def test_cancel(self):
        inv = self._obj_invoke_workflow()
        inv.cancel()
        assert inv.state in ["cancelled", "cancelling"]

    def test_wait(self):
        inv = self._obj_invoke_workflow()
        inv.wait()
        assert inv.state == "scheduled"

    def test_refresh(self):
        inv = self._obj_invoke_workflow()
        inv.state = "placeholder"
        # use wait_for_invocation() directly, because inv.wait() will update inv automatically
        self.gi.gi.invocations.wait_for_invocation(inv.id)
        inv.refresh()
        assert inv.state == "scheduled"

    def test_run_step_actions(self):
        inv = self.workflow_pause.invoke(
            inputs={"0": self.dataset},
            history=self.history,
        )
        for _ in range(20):
            with pytest.raises(bioblend.TimeoutException):
                inv.wait(maxwait=0.5, interval=0.5)
            inv.refresh()
            if len(inv.steps) >= 3:
                break
        assert inv.steps[2].action is None
        inv.run_step_actions([inv.steps[2]], [True])
        assert inv.steps[2].action is True

    def test_summary(self):
        inv = self._obj_invoke_workflow()
        inv.wait()
        summary = inv.summary()
        assert summary["populated_state"] == "ok"

    def test_step_jobs_summary(self):
        inv = self._obj_invoke_workflow()
        inv.wait()
        step_jobs_summary = inv.step_jobs_summary()
        assert len(step_jobs_summary) == 1
        assert step_jobs_summary[0]["populated_state"] == "ok"

    def test_report(self):
        inv = self._obj_invoke_workflow()
        report = inv.report()
        assert "paste_columns" in report["markdown"]

    @test_util.skip_unless_galaxy("release_20.09")
    def test_biocompute_object(self):
        inv = self._obj_invoke_workflow()
        inv.wait()
        biocompute_object = inv.biocompute_object()
        assert len(biocompute_object["description_domain"]["pipeline_steps"]) == 1

    def _obj_invoke_workflow(self) -> wrappers.Invocation:
        return self.workflow.invoke(
            inputs={"Input 1": self.dataset, "Input 2": self.dataset},
            history=self.history,
            inputs_by="name",
        )


@test_util.skip_unless_galaxy("release_19.09")
class TestObjInvocationClient(GalaxyObjectsTestBase):
    history: wrappers.History
    inv: wrappers.Invocation
    workflow: wrappers.Workflow

    @classmethod
    def setUpClass(cls) -> None:
        super().setUpClass()
        with open(SAMPLE_FN) as f:
            cls.workflow = cls.gi.workflows.import_new(f.read())
        cls.history = cls.gi.histories.create(name="TestGalaxyObjInvocationClient")
        dataset = cls.history.paste_content("1\t2\t3")
        cls.inv = cls.workflow.invoke(
            inputs={"Input 1": dataset, "Input 2": dataset},
            history=cls.history,
            inputs_by="name",
        )
        cls.inv.wait()

    @classmethod
    def tearDownClass(cls):
        cls.history.delete(purge=True)

    def test_get(self):
        inv = self.gi.invocations.get(self.inv.id)
        assert inv.id == self.inv.id
        assert inv.workflow_id == self.workflow.id
        assert inv.history_id == self.history.id
        assert inv.state == "scheduled"
        assert inv.update_time == self.inv.update_time
        assert inv.uuid == self.inv.uuid

    def test_get_previews(self):
        previews = self.gi.invocations.get_previews()
        assert {type(preview) for preview in previews} == {wrappers.InvocationPreview}
        inv_preview = next(p for p in previews if p.id == self.inv.id)
        assert inv_preview.id == self.inv.id
        assert inv_preview.workflow_id == self.workflow.id
        assert inv_preview.history_id == self.history.id
        assert inv_preview.state == "scheduled"
        assert inv_preview.update_time == self.inv.update_time
        assert inv_preview.uuid == self.inv.uuid

    def test_list(self):
        invs = self.gi.invocations.list()
        inv = next(i for i in invs if i.id == self.inv.id)
        assert inv.id == self.inv.id
        assert inv.workflow_id == self.workflow.id
        assert inv.history_id == self.history.id
        assert inv.state == "scheduled"
        assert inv.update_time == self.inv.update_time
        assert inv.uuid == self.inv.uuid
        assert len(self.inv.steps) > 0
        history = self.gi.histories.create(name="TestGalaxyObjInvocationClientList")
        assert self.gi.invocations.list(history=history) == []
        history.delete(purge=True)


class TestGalaxyInstance(GalaxyObjectsTestBase):
    def test_library(self):
        name = f"test_{uuid.uuid4().hex}"
        description, synopsis = "D", "S"
        lib = self.gi.libraries.create(name, description=description, synopsis=synopsis)
        assert lib.name == name
        assert lib.description == description
        assert lib.synopsis == synopsis
        assert len(lib.content_infos) == 1  # root folder
        assert len(lib.folder_ids) == 1
        assert len(lib.dataset_ids) == 0
        assert lib.id in [_.id for _ in self.gi.libraries.list()]
        lib.delete()
        assert not lib.is_mapped

    def test_workflow_from_str(self):
        with open(SAMPLE_FN) as f:
            wf = self.gi.workflows.import_new(f.read())
        self._check_and_del_workflow(wf)

    def test_workflow_collections_from_str(self):
        with open(SAMPLE_WF_COLL_FN) as f:
            wf = self.gi.workflows.import_new(f.read())
        self._check_and_del_workflow(wf)

    def test_workflow_parameter_input(self):
        with open(SAMPLE_WF_PARAMETER_INPUT_FN) as f:
            self.gi.workflows.import_new(f.read())

    def test_workflow_from_dict(self):
        with open(SAMPLE_FN) as f:
            wf = self.gi.workflows.import_new(json.load(f))
        self._check_and_del_workflow(wf)

    def test_workflow_publish_from_dict(self):
        with open(SAMPLE_FN) as f:
            wf = self.gi.workflows.import_new(json.load(f), publish=True)
        self._check_and_del_workflow(wf, check_is_public=True)

    def test_workflow_missing_tools(self):
        with open(SAMPLE_FN) as f:
            wf_dump = json.load(f)
        wf_info = self.gi.gi.workflows.import_workflow_dict(wf_dump)
        wf_dict = self.gi.gi.workflows.show_workflow(wf_info["id"])
        for id_, step in wf_dict["steps"].items():
            if step["type"] == "tool":
                for k in "tool_inputs", "tool_version":
                    wf_dict["steps"][id_][k] = None
        wf = wrappers.Workflow(wf_dict, gi=self.gi)
        assert not wf.is_runnable
        with pytest.raises(RuntimeError):
            wf.invoke()
        wf.delete()

    def test_workflow_export(self):
        with open(SAMPLE_FN) as f:
            wf1 = self.gi.workflows.import_new(f.read())
        wf2 = self.gi.workflows.import_new(wf1.export())
        assert wf1.id != wf2.id
        for wf in wf1, wf2:
            self._check_and_del_workflow(wf)

    def _check_and_del_workflow(self, wf: wrappers.Workflow, check_is_public: bool = False) -> None:
        # Galaxy appends additional text to imported workflow names
        assert wf.name.startswith("paste_columns")
        assert len(wf.steps) == 3
        for step_id, step in wf.steps.items():
            assert isinstance(step, wrappers.Step)
            assert step_id == step.id
            assert isinstance(step.tool_inputs, dict)
            if step.type == "tool":
                assert step.tool_id is not None
                assert step.tool_version is not None
                assert isinstance(step.input_steps, dict)
            elif step.type in ("data_collection_input", "data_input"):
                assert step.tool_id is None
                assert step.tool_version is None
                assert step.input_steps == {}
        wf_ids = {_.id for _ in self.gi.workflows.list()}
        assert wf.id in wf_ids
        if check_is_public:
            assert wf.published
        wf.delete()

    def test_workflow_from_shared(self):
        with open(SAMPLE_FN) as f:
            imported_wf = self.gi.workflows.import_new(json.load(f), publish=True)
        _, new_gi = test_util.new_user_gi(self.gi.gi)
        new_obj_gi = galaxy_instance.GalaxyInstance(
            new_gi.base_url,
            api_key=new_gi.key,
        )
        imported_wf_by_new_user = new_obj_gi.workflows.import_shared(imported_wf.id)
        assert isinstance(imported_wf_by_new_user, wrappers.Workflow)
        # The new StoredWorkflow's name is prepended with 'imported: ', but in
        # recent Galaxy versions the API returns the Workflow's name, which is not changed.
        assert imported_wf_by_new_user.name in {f"imported: {imported_wf.name}", imported_wf.name}
        assert not imported_wf_by_new_user.deleted
        assert not imported_wf_by_new_user.published
        imported_wf_by_new_user.delete()
        imported_wf.delete()

    def test_get_libraries(self):
        self._test_multi_get("libraries")

    def test_get_histories(self):
        self._test_multi_get("histories")

    def test_get_workflows(self):
        self._test_multi_get("workflows")

    def _normalized_functions(
        self, obj_type: Literal["histories", "libraries", "workflows"]
    ) -> tuple[Callable, dict[str, Any]]:
        if obj_type == "libraries":
            create: Callable = self.gi.libraries.create
            del_kwargs = {}
        elif obj_type == "histories":
            create = self.gi.histories.create
            del_kwargs = {"purge": True}
        elif obj_type == "workflows":

            def create(name):
                with open(SAMPLE_FN) as f:
                    d = json.load(f)
                d["name"] = name
                return self.gi.workflows.import_new(d)

            del_kwargs = {}
        return create, del_kwargs

    def _test_multi_get(self, obj_type: Literal["histories", "libraries", "workflows"]) -> None:
        obj_gi_client = getattr(self.gi, obj_type)
        create, del_kwargs = self._normalized_functions(obj_type)

        def ids(seq: Iterable[wrappers.Wrapper]) -> set[str]:
            return {_.id for _ in seq}

        names = [f"test_{uuid.uuid4().hex}" for _ in range(2)]
        objs = []
        try:
            objs = [create(_) for _ in names]
            assert ids(objs) <= ids(obj_gi_client.list())
            if obj_type != "workflows":
                filtered = obj_gi_client.list(name=names[0])
                assert len(filtered) == 1
                assert filtered[0].id == objs[0].id
                del_id = objs[-1].id
                objs.pop().delete(**del_kwargs)
                assert del_id in ids(obj_gi_client.get_previews(deleted=True))
            else:
                # Galaxy appends info strings to imported workflow names
                prev = obj_gi_client.get_previews()[0]
                filtered = obj_gi_client.list(name=prev.name)
                assert len(filtered) == 1
                assert filtered[0].id == prev.id
        finally:
            for o in objs:
                o.delete(**del_kwargs)

    def test_delete_libraries_by_name(self):
        self._test_delete_by_name("libraries")
        self._test_delete_by_ambiguous_name("libraries")

    def test_delete_histories_by_name(self):
        self._test_delete_by_name("histories")
        self._test_delete_by_ambiguous_name("histories")

    def test_delete_workflows_by_name(self):
        self._test_delete_by_name("workflows")
        self._test_delete_by_ambiguous_name("workflows")

    def _test_delete_by_name(self, obj_type: Literal["histories", "libraries", "workflows"]) -> None:
        obj_gi_client = getattr(self.gi, obj_type)
        create, del_kwargs = self._normalized_functions(obj_type)
        name = f"test_{uuid.uuid4().hex}"
        create(name)
        prevs = [_ for _ in obj_gi_client.get_previews(name=name) if not _.deleted]
        assert len(prevs) == 1
        del_kwargs["name"] = name
        obj_gi_client.delete(**del_kwargs)
        prevs = [_ for _ in obj_gi_client.get_previews(name=name) if not _.deleted]
        assert len(prevs) == 0

    def _test_delete_by_ambiguous_name(self, obj_type: Literal["histories", "libraries", "workflows"]) -> None:
        obj_gi_client = getattr(self.gi, obj_type)
        create, del_kwargs = self._normalized_functions(obj_type)
        name = f"test_{uuid.uuid4().hex}"
        objs = [create(name) for _ in range(2)]
        prevs = [_ for _ in obj_gi_client.get_previews(name=name) if not _.deleted]
        assert len(prevs) == len(objs)
        del_kwargs["name"] = name
        with pytest.raises(ValueError):
            obj_gi_client.delete(**del_kwargs)
        # Cleanup
        del del_kwargs["name"]
        for prev in prevs:
            del_kwargs["id_"] = prev.id
            obj_gi_client.delete(**del_kwargs)


class TestLibrary(GalaxyObjectsTestBase):
    # just something that can be expected to be always up
    DS_URL = "https://tools.ietf.org/rfc/rfc1866.txt"

    def setUp(self):
        super().setUp()
        self.lib = self.gi.libraries.create(f"test_{uuid.uuid4().hex}")

    def tearDown(self):
        self.lib.delete()

    def test_root_folder(self):
        r = self.lib.root_folder
        assert r.parent is None

    def test_folder(self):
        name = f"test_{uuid.uuid4().hex}"
        desc = "D"
        folder = self.lib.create_folder(name, description=desc)
        assert folder.name == name
        assert folder.description == desc
        assert folder.container is self.lib
        assert folder.parent is not None
        assert folder.parent.id == self.lib.root_folder.id
        assert len(self.lib.content_infos) == 2
        assert len(self.lib.folder_ids) == 2
        assert folder.id in self.lib.folder_ids
        retrieved = self.lib.get_folder(folder.id)
        assert folder.id == retrieved.id

    def _check_datasets(self, dss: Collection[wrappers.LibraryDataset]) -> None:
        assert len(dss) == len(self.lib.dataset_ids)
        assert {_.id for _ in dss} == set(self.lib.dataset_ids)
        for ds in dss:
            assert ds.container is self.lib

    def test_dataset(self):
        folder = self.lib.create_folder(f"test_{uuid.uuid4().hex}")
        ds = self.lib.upload_data(FOO_DATA, folder=folder)
        assert len(self.lib.content_infos) == 3
        assert len(self.lib.folder_ids) == 2
        self._check_datasets([ds])

    def test_dataset_from_url(self):
        if is_reachable(self.DS_URL):
            ds = self.lib.upload_from_url(self.DS_URL)
            self._check_datasets([ds])
        else:
            self.skipTest(f"{self.DS_URL} not reachable")

    def test_dataset_from_local(self):
        with tempfile.NamedTemporaryFile(mode="w", prefix="bioblend_test_") as f:
            f.write(FOO_DATA)
            f.flush()
            ds = self.lib.upload_from_local(f.name)
        self._check_datasets([ds])

    def test_datasets_from_fs(self):
        bnames = [f"f{i}.txt" for i in range(2)]
        dss, fnames = upload_from_fs(self.lib, bnames)
        self._check_datasets(dss)
        dss, fnames = upload_from_fs(self.lib, bnames, link_data_only="link_to_files")
        for ds, fn in zip(dss, fnames):
            assert ds.file_name == fn

    def test_copy_from_dataset(self):
        hist = self.gi.histories.create(f"test_{uuid.uuid4().hex}")
        try:
            hda = hist.paste_content(FOO_DATA)
            ds = self.lib.copy_from_dataset(hda)
        finally:
            hist.delete(purge=True)
        self._check_datasets([ds])

    def test_get_dataset(self):
        ds = self.lib.upload_data(FOO_DATA)
        retrieved = self.lib.get_dataset(ds.id)
        assert ds.id == retrieved.id

    def test_get_datasets(self):
        bnames = [f"f{i}.txt" for i in range(2)]
        dss, _ = upload_from_fs(self.lib, bnames)
        retrieved = self.lib.get_datasets()
        assert len(dss) == len(retrieved)
        assert {_.id for _ in dss} == {_.id for _ in retrieved}
        name = f"/{bnames[0]}"
        selected = self.lib.get_datasets(name=name)
        assert len(selected) == 1
        assert selected[0].name == bnames[0]


class TestLDContents(GalaxyObjectsTestBase):
    def setUp(self):
        super().setUp()
        self.lib = self.gi.libraries.create(f"test_{uuid.uuid4().hex}")
        self.ds = self.lib.upload_data(FOO_DATA)
        self.ds.wait()

    def tearDown(self):
        self.lib.delete()

    def test_dataset_get_stream(self):
        for idx, c in enumerate(self.ds.get_stream(chunk_size=1)):
            assert FOO_DATA[idx].encode() == c

    def test_dataset_peek(self):
        fetched_data = self.ds.peek(chunk_size=4)
        assert FOO_DATA[0:4].encode() == fetched_data

    def test_dataset_download(self):
        with tempfile.TemporaryFile() as f:
            self.ds.download(f)
            f.seek(0)
            assert FOO_DATA.encode() == f.read()

    def test_dataset_get_contents(self):
        assert FOO_DATA.encode() == self.ds.get_contents()

    def test_dataset_delete(self):
        self.ds.delete()
        # Cannot test this yet because the 'deleted' attribute is not exported
        # by the API at the moment
        # assert self.ds.deleted

    def test_dataset_update(self):
        new_name = f"test_{uuid.uuid4().hex}"
        new_misc_info = f"Annotation for {new_name}"
        new_genome_build = "hg19"
        updated_ldda = self.ds.update(name=new_name, misc_info=new_misc_info, genome_build=new_genome_build)
        assert self.ds.id == updated_ldda.id
        assert self.ds.name == new_name
        assert self.ds.misc_info == new_misc_info
        assert self.ds.genome_build == new_genome_build


class TestHistory(GalaxyObjectsTestBase):
    def setUp(self):
        super().setUp()
        self.hist = self.gi.histories.create(f"test_{uuid.uuid4().hex}")

    def tearDown(self):
        self.hist.delete(purge=True)

    def test_create_delete(self):
        name = f"test_{uuid.uuid4().hex}"
        hist = self.gi.histories.create(name)
        assert hist.name == name
        hist_id = hist.id
        assert hist_id in [_.id for _ in self.gi.histories.list()]
        hist.delete(purge=True)
        assert not hist.is_mapped
        h = self.gi.histories.get(hist_id)
        assert h.deleted

    def _check_dataset(self, hda: wrappers.HistoryDatasetAssociation) -> None:
        assert isinstance(hda, wrappers.HistoryDatasetAssociation)
        assert hda.container is self.hist
        assert len(self.hist.dataset_ids) == 1
        assert self.hist.dataset_ids[0] == hda.id

    def test_import_dataset(self):
        lib = self.gi.libraries.create(f"test_{uuid.uuid4().hex}")
        lds = lib.upload_data(FOO_DATA)
        assert len(self.hist.dataset_ids) == 0
        hda = self.hist.import_dataset(lds)
        lib.delete()
        self._check_dataset(hda)

    def test_upload_file(self):
        with tempfile.NamedTemporaryFile(mode="w", prefix="bioblend_test_") as f:
            f.write(FOO_DATA)
            f.flush()
            hda = self.hist.upload_file(f.name)
        self._check_dataset(hda)

    def test_paste_content(self):
        hda = self.hist.paste_content(FOO_DATA)
        self._check_dataset(hda)

    def test_get_dataset(self):
        hda = self.hist.paste_content(FOO_DATA)
        retrieved = self.hist.get_dataset(hda.id)
        assert hda.id == retrieved.id

    def test_get_datasets(self):
        bnames = [f"f{i}.txt" for i in range(2)]
        lib = self.gi.libraries.create(f"test_{uuid.uuid4().hex}")
        lds = upload_from_fs(lib, bnames)[0]
        hdas = [self.hist.import_dataset(_) for _ in lds]
        lib.delete()
        retrieved = self.hist.get_datasets()
        assert len(hdas) == len(retrieved)
        assert {_.id for _ in hdas} == {_.id for _ in retrieved}
        selected = self.hist.get_datasets(name=bnames[0])
        assert len(selected) == 1
        assert selected[0].name == bnames[0]

    def test_export_and_download(self):
        jeha_id = self.hist.export(wait=True, maxwait=60)
        assert jeha_id
        tempdir = tempfile.mkdtemp(prefix="bioblend_test_")
        temp_fn = os.path.join(tempdir, "export.tar.gz")
        try:
            with open(temp_fn, "wb") as fo:
                self.hist.download(jeha_id, fo)
            assert tarfile.is_tarfile(temp_fn)
        finally:
            shutil.rmtree(tempdir)

    def test_update(self):
        new_name = f"test_{uuid.uuid4().hex}"
        new_annotation = f"Annotation for {new_name}"
        new_tags = ["tag1", "tag2"]
        updated_hist = self.hist.update(name=new_name, annotation=new_annotation, tags=new_tags)
        assert self.hist.id == updated_hist.id
        assert self.hist.name == new_name
        assert self.hist.annotation == new_annotation
        assert self.hist.tags == new_tags
        updated_hist = self.hist.update(published=True)
        assert self.hist.id == updated_hist.id
        assert self.hist.published

    def test_create_dataset_collection(self):
        self._create_collection_description()
        hdca = self.hist.create_dataset_collection(self.collection_description, copy_elements=False)
        assert isinstance(hdca, wrappers.HistoryDatasetCollectionAssociation)
        assert hdca.collection_type == "list"
        assert hdca.container is self.hist
        assert len(hdca.elements) == 2
        assert self.dataset1.id == hdca.elements[0]["object"]["id"]
        assert self.dataset2.id == hdca.elements[1]["object"]["id"]

    def test_delete_dataset_collection(self):
        self._create_collection_description()
        hdca = self.hist.create_dataset_collection(self.collection_description, copy_elements=False)
        hdca.delete()
        assert hdca.deleted

    def _create_collection_description(self) -> None:
        self.dataset1 = self.hist.paste_content(FOO_DATA)
        self.dataset2 = self.hist.paste_content(FOO_DATA_2)
        self.collection_description = dataset_collections.CollectionDescription(
            name="MyDatasetList",
            elements=[
                dataset_collections.HistoryDatasetElement(name="sample1", id=self.dataset1.id),
                dataset_collections.HistoryDatasetElement(name="sample2", id=self.dataset2.id),
            ],
        )


class TestHDAContents(GalaxyObjectsTestBase):
    def setUp(self):
        super().setUp()
        self.hist = self.gi.histories.create(f"test_{uuid.uuid4().hex}")
        self.ds = self.hist.paste_content(FOO_DATA)

    def tearDown(self):
        self.hist.delete(purge=True)

    def test_dataset_get_stream(self):
        for idx, c in enumerate(self.ds.get_stream(chunk_size=1)):
            assert FOO_DATA[idx].encode() == c

    def test_dataset_peek(self):
        fetched_data = self.ds.peek(chunk_size=4)
        assert FOO_DATA[0:4].encode() == fetched_data

    def test_dataset_download(self):
        with tempfile.TemporaryFile() as f:
            self.ds.download(f)
            f.seek(0)
            assert FOO_DATA.encode() == f.read()

    def test_dataset_get_contents(self):
        assert FOO_DATA.encode() == self.ds.get_contents()

    def test_dataset_update(self):
        new_name = f"test_{uuid.uuid4().hex}"
        new_annotation = f"Annotation for {new_name}"
        new_genome_build = "hg19"
        updated_hda = self.ds.update(name=new_name, annotation=new_annotation, genome_build=new_genome_build)
        assert self.ds.id == updated_hda.id
        assert self.ds.name == new_name
        assert self.ds.annotation == new_annotation
        assert self.ds.genome_build == new_genome_build

    def test_dataset_delete(self):
        self.ds.delete()
        assert self.ds.deleted
        assert not self.ds.purged

    def test_dataset_purge(self):
        self.ds.delete(purge=True, wait=True)
        assert self.ds.deleted
        assert self.ds.purged


@test_util.skip_unless_galaxy("release_19.09")
class TestRunWorkflow(GalaxyObjectsTestBase):
    def setUp(self):
        super().setUp()
        self.lib = self.gi.libraries.create(f"test_{uuid.uuid4().hex}")
        with open(SAMPLE_FN) as f:
            self.wf = self.gi.workflows.import_new(f.read())
        self.contents = ["one\ntwo\n", "1\n2\n"]
        self.inputs = [self.lib.upload_data(_) for _ in self.contents]

    def tearDown(self):
        self.wf.delete()
        self.lib.delete()

    def _test(self, existing_hist: bool = False, pass_params: bool = False) -> None:
        hist_name = f"test_{uuid.uuid4().hex}"
        if existing_hist:
            hist: Union[str, wrappers.History] = self.gi.histories.create(hist_name)
        else:
            hist = hist_name
        if pass_params:
            params = {"Paste1": {"delimiter": "U"}}
            sep = "_"  # 'U' maps to '_' in the paste tool
        else:
            params = None
            sep = "\t"  # default
        input_map = {"Input 1": self.inputs[0], "Input 2": self.inputs[1]}
        sys.stderr.write(os.linesep)
        inv = self.wf.invoke(inputs=input_map, params=params, history=hist, inputs_by="name")
        out_hist = self.gi.histories.get(inv.history_id)
        inv.wait()
        last_step = inv.sorted_steps_by()[-1]
        out_ds = last_step.get_outputs()["out_file1"]
        assert out_ds.container.id == out_hist.id
        res = out_ds.get_contents()
        exp_rows = zip(*(_.splitlines() for _ in self.contents))
        exp_res = ("\n".join(sep.join(t) for t in exp_rows) + "\n").encode()
        assert res == exp_res
        if isinstance(hist, wrappers.History):  # i.e. existing_hist == True
            assert out_hist.id == hist.id
        out_hist.delete(purge=True)

    def test_existing_history(self) -> None:
        self._test(existing_hist=True)

    def test_new_history(self) -> None:
        self._test(existing_hist=False)

    def test_params(self) -> None:
        self._test(pass_params=True)


@test_util.skip_unless_galaxy("release_19.09")
class TestRunDatasetCollectionWorkflow(GalaxyObjectsTestBase):
    def setUp(self):
        super().setUp()
        with open(SAMPLE_WF_COLL_FN) as f:
            self.wf = self.gi.workflows.import_new(f.read())
        self.hist = self.gi.histories.create(f"test_{uuid.uuid4().hex}")

    def tearDown(self):
        self.wf.delete()
        self.hist.delete(purge=True)

    def test_run_workflow_with_dataset_collection(self):
        dataset1 = self.hist.paste_content(FOO_DATA)
        dataset2 = self.hist.paste_content(FOO_DATA_2)
        collection_description = dataset_collections.CollectionDescription(
            name="MyDatasetList",
            elements=[
                dataset_collections.HistoryDatasetElement(name="sample1", id=dataset1.id),
                dataset_collections.HistoryDatasetElement(name="sample2", id=dataset2.id),
            ],
        )
        dataset_collection = self.hist.create_dataset_collection(collection_description, copy_elements=False)
        assert len(self.hist.content_infos) == 3
        input_map = {"0": dataset_collection, "1": dataset1}
        inv = self.wf.invoke(input_map, history=self.hist)
        inv.wait()
        self.hist.refresh()
        assert len(self.hist.content_infos) == 6
        last_step = inv.sorted_steps_by()[-1]
        out_hdca = last_step.get_output_collections()["out_file1"]
        assert out_hdca.collection_type == "list"
        assert len(out_hdca.elements) == 2
        assert out_hdca.container.id == self.hist.id


class TestJob(GalaxyObjectsTestBase):
    def test_get(self):
        job_prevs = self.gi.jobs.get_previews()
        if len(job_prevs) > 0:
            job_prev = job_prevs[0]
            assert isinstance(job_prev, wrappers.JobPreview)
            job = self.gi.jobs.get(job_prev.id)
            assert isinstance(job, wrappers.Job)
            assert job.id == job_prev.id
        for job in self.gi.jobs.list():
            assert isinstance(job, wrappers.Job)
