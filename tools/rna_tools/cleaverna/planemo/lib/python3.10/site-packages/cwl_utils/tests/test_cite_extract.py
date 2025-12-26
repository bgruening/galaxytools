# SPDX-License-Identifier: Apache-2.0
"""Tests for cwl-cite-extract."""
import pytest

from cwl_utils.cite_extract import arg_parser, run

from .util import get_data


def test_cite_extract_simple(capsys: pytest.CaptureFixture[str]) -> None:
    """Test the citation extraction, simply."""
    assert run(arg_parser().parse_args([get_data("testdata/seqtk_seq.cwl")])) == 0
    captured = capsys.readouterr()
    assert captured.out == "Package: seqtk, version: ['r93'], specs: None\n"
    assert captured.err == ""


def test_cite_extract_workflow_no_results(capsys: pytest.CaptureFixture[str]) -> None:
    """Attempt to extract citations from a workflow without any SoftwareRequirements."""
    assert (
        run(
            arg_parser().parse_args([get_data("testdata/checker_wf/functional-wf.cwl")])
        )
        == 0
    )
    captured = capsys.readouterr()
    assert captured.out == ""
    assert captured.err == ""
