"""Tests related to the --dir command line option."""

import filecmp
from pathlib import Path

from cwlupgrader.main import main

from .util import get_data


def test_draft3_workflow(tmp_path: Path) -> None:
    """Confirm that --dir works when the directory doesn't exist yet."""
    out_dir = tmp_path / "new"
    main([f"--dir={out_dir}", "--v1-only", get_data("testdata/draft-3/wf.cwl")])
    result = filecmp.cmp(
        get_data("testdata/v1.0/wf.cwl"),
        out_dir / "wf.cwl",
        shallow=False,
    )
    assert result
