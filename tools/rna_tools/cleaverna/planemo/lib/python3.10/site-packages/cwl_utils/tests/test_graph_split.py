# SPDX-License-Identifier: Apache-2.0
"""Test the CWL $graph document splitter tool."""
import json
from io import StringIO
from pathlib import Path

import pytest
import requests
from cwltool.tests.util import get_main_output

from cwl_utils.graph_split import graph_split

from .util import get_path

URI = (
    "https://gist.githubusercontent.com/altairwei/"
    "6a0097db95cad23de36f825ed3b9f4b0/raw/"
    "83f332931c3093ee73554cd7f60054ce17d03239/rhapsody_wta_1.8.packed.cwl"
)


def test_graph_split(tmp_path: Path) -> None:
    """Confirm that a user provided example produces no exception."""
    sourceIO = StringIO(requests.get(URI).text)
    sourceIO.name = URI
    graph_split(sourceIO, tmp_path, "yaml", "main.cwl", True)


def test_graph_split_offline(tmp_path: Path) -> None:
    """Confirm that a local provided example produces no exception."""
    with get_path("testdata/js-expr-req-wf.cwl").open() as handle:
        graph_split(handle, tmp_path, "yaml", "main.cwl", True)
    target = tmp_path / "wf.cwl"
    assert target.exists()
    code, stdout, stderr = get_main_output(["--debug", str(target)])
    assert code == 0, stderr
    assert (
        json.loads(stdout)["out"]["checksum"]
        == "sha1$7448d8798a4380162d4b56f9b452e2f6f9e24e7a"
    )


def test_graph_split_json_offline(tmp_path: Path) -> None:
    """Confirm that a local provided example produces no exception in JSON mode."""
    target = tmp_path / "subdir" / "wf.cwl"
    with get_path("testdata/js-expr-req-wf.cwl").open() as handle:
        graph_split(handle, target.parent, "json", "main.cwl", True)
    assert target.exists()
    code, stdout, stderr = get_main_output(["--debug", str(target)])
    assert code == 0, stderr
    assert (
        json.loads(stdout)["out"]["checksum"]
        == "sha1$7448d8798a4380162d4b56f9b452e2f6f9e24e7a"
    )


def test_graph_split_bad_path() -> None:
    """Expect an exception when the target directory parent does not exist."""
    with get_path("testdata/js-expr-req-wf.cwl").open() as handle:
        with pytest.raises(NotADirectoryError):
            graph_split(
                handle, Path("/__non_existent/tmp_path"), "json", "main.cwl", True
            )


def test_graph_split_complex1(tmp_path: Path) -> None:
    """Split a more complex graph with SchemaDefRequirement and $import."""
    with get_path("testdata/remote-cwl/wf1-packed.cwl").open() as handle:
        graph_split(handle, tmp_path, "yaml", "main.cwl", False)


def test_graph_split_complex2(tmp_path: Path) -> None:
    """Split another complex graph with SchemaDefRequirement and $import."""
    with get_path("testdata/workflows/wf5-packed.cwl").open() as handle:
        graph_split(handle, tmp_path, "yaml", "main.cwl", False)
