import filecmp
from pathlib import Path

from cwlupgrader.main import load_cwl_document, main, upgrade_document

from .util import get_data


def test_draft3_workflow(tmp_path: Path) -> None:
    """Basic draft3 to CWL v1.1 test."""
    main([f"--dir={tmp_path}", "--v1-only", get_data("testdata/draft-3/wf.cwl")])
    result = filecmp.cmp(
        get_data("testdata/v1.0/wf.cwl"),
        tmp_path / "wf.cwl",
        shallow=False,
    )
    assert result


def test_draft3_tool_long_form_arrays(tmp_path: Path) -> None:
    """Draft-3 document with long form array inputs."""
    main(
        [
            f"--dir={tmp_path}",
            "--v1-only",
            get_data("testdata/draft-3/attributor-prok-cheetah.cwl"),
        ]
    )
    result = filecmp.cmp(
        get_data("testdata/v1.0/attributor-prok-cheetah.cwl"),
        tmp_path / "attributor-prok-cheetah.cwl",
        shallow=False,
    )
    assert result


def test_invalid_target(tmp_path: Path) -> None:
    """Test for invalid target version"""
    doc = load_cwl_document(get_data("testdata/v1.0/listing_deep1.cwl"))
    result = upgrade_document(doc, str(tmp_path), "invalid-version")
    assert result is None


def test_v1_0_to_v1_1_load_listing(tmp_path: Path) -> None:
    """Basic CWL v1.0 to CWL v1.1 test with LoadListingRequirement (map notation)."""
    doc = load_cwl_document(get_data("testdata/v1.0/listing_deep1.cwl"))
    upgraded = upgrade_document(doc, str(tmp_path), "v1.1")
    expected = load_cwl_document(get_data("testdata/v1.1/listing_deep1.cwl"))
    assert upgraded == expected


def test_v1_0_to_v1_1_load_listing_arr(tmp_path: Path) -> None:
    """Basic CWL v1.0 to CWL v1.1 test with LoadListingRequirement (array notation)."""
    doc = load_cwl_document(get_data("testdata/v1.0/listing_deep1-arr.cwl"))
    upgraded = upgrade_document(doc, str(tmp_path), "v1.1")
    expected = load_cwl_document(get_data("testdata/v1.1/listing_deep1-arr.cwl"))
    assert upgraded == expected


def test_v1_0_to_v1_1_network_access(tmp_path: Path) -> None:
    """Basic CWL v1.0 to CWL v1.1 test with NetworkAccess."""
    doc = load_cwl_document(get_data("testdata/v1.0/networkaccess.cwl"))
    upgraded = upgrade_document(doc, str(tmp_path), "v1.1")
    expected = load_cwl_document(get_data("testdata/v1.1/networkaccess.cwl"))
    assert upgraded == expected


def test_v1_1_to_v1_2(tmp_path: Path) -> None:
    """Basic CWL v1.1 to CWL v1.2 test."""
    doc = load_cwl_document(get_data("testdata/v1.1/listing_deep1.cwl"))
    upgraded = upgrade_document(doc, str(tmp_path), "v1.2")
    expected = load_cwl_document(get_data("testdata/v1.2/listing_deep1.cwl"))
    assert upgraded == expected


def test_v1_2_to_v1_2(tmp_path: Path) -> None:
    """CWL v1.2 to CWL v1.2 no change test."""
    doc = load_cwl_document(get_data("testdata/v1.2/networkaccess.cwl"))
    upgraded = upgrade_document(doc, str(tmp_path), "v1.2")
    expected = load_cwl_document(get_data("testdata/v1.2/networkaccess.cwl"))
    assert upgraded == expected


def test_v1_2_to_latest(tmp_path: Path) -> None:
    """CWL v1.2 to latest no change test."""
    doc = load_cwl_document(get_data("testdata/v1.2/networkaccess.cwl"))
    upgraded = upgrade_document(doc, str(tmp_path), "latest")
    expected = load_cwl_document(get_data("testdata/v1.2/networkaccess.cwl"))
    assert upgraded == expected


def test_packed_graph(tmp_path: Path) -> None:
    """Test packed document with $graph."""
    main(
        [f"--dir={tmp_path}", "--v1.1-only", get_data("testdata/v1.0/conflict-wf.cwl")]
    )
    assert filecmp.cmp(
        get_data("testdata/v1.1/conflict-wf.cwl"),
        tmp_path / "conflict-wf.cwl",
        shallow=False,
    )


def test_multi_version_upgrade_external_steps(tmp_path: Path) -> None:
    """Test 1.0 to 1.2 upgrade of Workflow with external steps."""
    main([f"--dir={tmp_path}", get_data("testdata/v1.0/1st-workflow.cwl")])
    assert filecmp.cmp(
        get_data("testdata/v1.2/arguments.cwl"),
        tmp_path / "arguments.cwl",
        shallow=False,
    )
    assert filecmp.cmp(
        get_data("testdata/v1.2/tar-param.cwl"),
        tmp_path / "tar-param.cwl",
        shallow=False,
    )
