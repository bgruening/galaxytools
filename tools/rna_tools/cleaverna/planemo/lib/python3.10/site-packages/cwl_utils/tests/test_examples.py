# SPDX-License-Identifier: Apache-2.0
"""Tests of example Python scripts."""
import os
import runpy
from pathlib import Path


def test_load_example() -> None:
    """Test the load_cwl_by_path.py script."""
    cwd = Path.cwd()
    parent = Path(__file__).resolve().parent
    os.chdir(parent.parent)
    result_raw = runpy.run_path(str(parent / "load_cwl_by_path.py"))
    os.chdir(cwd)
    result = result_raw["saved_obj"]
    assert result["class"] == "Workflow"
