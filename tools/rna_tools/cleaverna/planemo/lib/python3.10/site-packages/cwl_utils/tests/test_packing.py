from typing import Any

from cwl_utils.pack import pack

from .util import get_data


def _find(l_item: list[Any], key: str, val: str) -> Any:
    return next(_x for _x in l_item if _x[key] == val)


def test_port_normalization() -> None:
    cwl = pack(get_data("testdata/remote-cwl/wf1.cwl"))
    step_s1 = _find(cwl.get("steps", []), "id", "s1")
    step_in1 = _find(step_s1.get("in"), "id", "in1")
    assert step_in1["source"] == "in1"

    cwl = pack(get_data("testdata/wf2.cwl"))
    step_s1 = _find(cwl.get("steps", []), "id", "s1")
    step_in1 = _find(step_s1.get("in"), "id", "in1")
    assert step_in1["source"] == "in1"

    out1 = _find(cwl.get("outputs", []), "id", "out1")
    assert out1.get("outputSource") == "s2/out1"


def test_include() -> None:
    cwl = pack(get_data("testdata/remote-cwl/tool1.cwl"))
    assert "arguments" in cwl
    assert isinstance(cwl.get("arguments"), list)

    inline_js_req = _find(
        cwl.get("requirements", []), "class", "InlineJavascriptRequirement"
    )
    include_js = inline_js_req.get("expressionLib")

    assert isinstance(include_js, list)
    assert "engineers walk into a" in include_js[0]


def test_schema_def1() -> None:
    cwl = pack(get_data("testdata/remote-cwl/tool2.cwl"))
    _type = _find(cwl.get("inputs", []), "id", "in1").get("type")
    assert isinstance(_type, dict)
    assert _type.get("type") == "array"


def test_schema_def2() -> None:
    cwl = pack(get_data("testdata/wf2.cwl"))
    _type = _find(cwl.get("inputs", []), "id", "in2").get("type")
    assert isinstance(_type, dict)
    assert _type.get("type") == "enum"


def test_step_packing() -> None:
    cwl = pack(get_data("testdata/remote-cwl/wf1.cwl"))
    s1 = _find(cwl.get("steps", []), "id", "s1")
    tool2 = s1.get("run")
    _type = _find(tool2.get("inputs"), "id", "in1").get("type")
    assert isinstance(_type, dict)
    assert _type.get("type") == "array"


def test_embedded_packing() -> None:
    pack(get_data("testdata/workflows/count-lines16-wf.cwl"))


def test_remote_packing() -> None:
    cwl = pack(
        "https://raw.githubusercontent.com/kaushik-work/sbpack/master/tests/wf2.cwl"
    )
    s1 = _find(cwl.get("steps", []), "id", "s1")
    wf1 = s1.get("run")
    assert wf1.get("class") == "Workflow"

    tool2 = _find(wf1.get("steps"), "id", "s1").get("run")
    _type = _find(tool2.get("inputs"), "id", "in1").get("type")
    assert isinstance(_type, dict)
    assert _type.get("type") == "array"


def test_remote_packing_github_soft_links() -> None:
    cwl = pack(
        "https://raw.githubusercontent.com/rabix/sbpack/master/tests/workflows/wf5.cwl"
    )
    s1 = _find(cwl.get("steps", []), "id", "s1")
    tool1 = s1.get("run")
    assert tool1.get("class") == "CommandLineTool"


def test_already_packed_graph() -> None:
    """Workflow already packed in a $graph."""
    cwl = pack(get_data("testdata/workflows/scatter-wf4.cwl"))
    assert "inputs" not in cwl
    assert "outputs" not in cwl
    assert "$graph" in cwl
    assert "requirements" not in cwl


def test_import_in_type() -> None:
    cwl = pack(get_data("testdata/workflows/import-in-type.cwl"))
    assert cwl["inputs"][0]["type"] == ["File", "Directory"]
