# SPDX-License-Identifier: Apache-2.0
"""Test the CWL Expression refactoring tool."""
import os
import shutil
import sys
import tarfile
from collections.abc import Generator
from pathlib import Path
from typing import TYPE_CHECKING, cast

import pytest
import requests
from _pytest.tmpdir import TempPathFactory
from pytest import raises

import cwl_utils.parser.cwl_v1_0 as parser
import cwl_utils.parser.cwl_v1_1 as parser1
import cwl_utils.parser.cwl_v1_2 as parser2
from cwl_utils.cwl_v1_0_expression_refactor import traverse as traverse0
from cwl_utils.cwl_v1_1_expression_refactor import traverse as traverse1
from cwl_utils.cwl_v1_2_expression_refactor import traverse as traverse2
from cwl_utils.errors import WorkflowException
from cwl_utils.expression_refactor import run as expression_refactor

from .util import get_data

if TYPE_CHECKING:
    from http.client import HTTPResponse


def test_v1_0_workflow_top_level_format_expr() -> None:
    """Test for the correct error when converting a format expression in a workflow level input."""
    with raises(WorkflowException, match=r".*format specification.*"):
        result, modified = traverse0(
            parser.load_document(get_data("testdata/workflow_input_format_expr.cwl")),
            False,
            False,
            False,
            False,
        )


def test_v1_0_workflow_top_level_sf_expr() -> None:
    """Test for the correct error when converting a secondaryFiles expression in a workflow level input."""
    with raises(WorkflowException, match=r".*secondaryFiles.*"):
        result, modified = traverse0(
            parser.load_document(get_data("testdata/workflow_input_sf_expr.cwl")),
            False,
            False,
            False,
            False,
        )


def test_v1_0_workflow_top_level_sf_expr_array() -> None:
    """Test correct error when converting a secondaryFiles expression (array form) in a workflow level input."""  # noqa: B950
    with raises(WorkflowException, match=r".*secondaryFiles.*"):
        result, modified = traverse0(
            parser.load_document(get_data("testdata/workflow_input_sf_expr_array.cwl")),
            False,
            False,
            False,
            False,
        )


def test_v1_1_workflow_top_level_format_expr() -> None:
    """Test for the correct error when converting a format expression in a workflow level input."""
    with raises(WorkflowException, match=r".*format specification.*"):
        result, modified = traverse1(
            parser1.load_document(
                get_data("testdata/workflow_input_format_expr_v1_1.cwl")
            ),
            False,
            False,
            False,
            False,
        )


def test_v1_1_workflow_top_level_sf_expr() -> None:
    """Test for the correct error when converting a secondaryFiles expression in a workflow level input."""
    with raises(WorkflowException, match=r".*secondaryFiles.*"):
        result, modified = traverse1(
            parser1.load_document(get_data("testdata/workflow_input_sf_expr_v1_1.cwl")),
            False,
            False,
            False,
            False,
        )


def test_v1_1_workflow_top_level_sf_expr_array() -> None:
    """Test for the correct error when converting a secondaryFiles expression (array form) in a workflow level input."""  # noqa: B950
    with raises(WorkflowException, match=r".*secondaryFiles.*"):
        result, modified = traverse1(
            parser1.load_document(
                get_data("testdata/workflow_input_sf_expr_array_v1_1.cwl")
            ),
            False,
            False,
            False,
            False,
        )


def test_v1_2_workflow_top_level_format_expr() -> None:
    """Test for the correct error when converting a format expression in a workflow level input."""
    with raises(WorkflowException, match=r".*format specification.*"):
        result, modified = traverse2(
            parser2.load_document(
                get_data("testdata/workflow_input_format_expr_v1_2.cwl")
            ),
            False,
            False,
            False,
            False,
        )


def test_v1_2_workflow_top_level_sf_expr() -> None:
    """Test for the correct error when converting a secondaryFiles expression in a workflow level input."""
    with raises(WorkflowException, match=r".*secondaryFiles.*"):
        result, modified = traverse2(
            parser2.load_document(get_data("testdata/workflow_input_sf_expr_v1_2.cwl")),
            False,
            False,
            False,
            False,
        )


def test_v1_2_workflow_top_level_sf_expr_array() -> None:
    """Test for the correct error when converting a secondaryFiles expression (array form) in a workflow level input."""  # noqa: B950
    with raises(WorkflowException, match=r".*secondaryFiles.*"):
        result, modified = traverse2(
            parser2.load_document(
                get_data("testdata/workflow_input_sf_expr_array_v1_2.cwl")
            ),
            False,
            False,
            False,
            False,
        )


def test_v1_0_step_valuefrom_expr_multisource() -> None:
    """Convert a valueFrom expression that has multiple sources."""
    result, modified = traverse0(
        parser.load_document(get_data("testdata/step-valuefrom2-wf_v1_0.cwl")),
        False,
        False,
        False,
        False,
    )


def test_v1_1_step_valuefrom_expr_multisource() -> None:
    """Convert a valueFrom expression that has multiple sources."""
    result, modified = traverse1(
        parser1.load_document(get_data("testdata/step-valuefrom2-wf_v1_1.cwl")),
        False,
        False,
        False,
        False,
    )


def test_v1_2_step_valuefrom_expr_multisource() -> None:
    """Convert a valueFrom expression that has multiple sources."""
    result, modified = traverse2(
        parser2.load_document(get_data("testdata/step-valuefrom2-wf_v1_2.cwl")),
        False,
        False,
        False,
        False,
    )


def test_v1_0_step_valuefrom_expr_sibling_inputs() -> None:
    """Convert a valueFrom expression from a step input that has uninvolved sibling inputs."""
    result, modified = traverse0(
        parser.load_document(get_data("testdata/step-valuefrom3-wf_v1_0.cwl")),
        False,
        False,
        False,
        False,
    )


def test_v1_1_step_valuefrom_expr_sibling_inputs() -> None:
    """Convert a valueFrom expression from a step input that has uninvolved sibling inputs."""
    result, modified = traverse1(
        parser1.load_document(get_data("testdata/step-valuefrom3-wf_v1_1.cwl")),
        False,
        False,
        False,
        False,
    )


def test_v1_2_step_valuefrom_expr_sibling_inputs() -> None:
    """Convert a valueFrom expression from a step input that has uninvolved sibling inputs."""
    result, modified = traverse2(
        parser2.load_document(get_data("testdata/step-valuefrom3-wf_v1_2.cwl")),
        False,
        False,
        False,
        False,
    )


def test_v1_2_workflow_output_pickvalue_expr() -> None:
    """Convert a workflow output pickValue expression."""
    result, modified = traverse2(
        parser2.load_document(get_data("testdata/cond-wf-003.1.cwl")),
        False,
        False,
        False,
        False,
    )


def test_expression_refactor(tmp_path: Path) -> None:
    """Functional test."""
    input_path = get_data("testdata/cond-wf-003.1.cwl")
    result = expression_refactor([str(tmp_path), input_path])
    assert result == 0


def test_expression_refactor_noop_solo(tmp_path: Path) -> None:
    """Functional test."""
    input_path = get_data("testdata/dockstore-tool-md5sum.cwl")
    result = expression_refactor([str(tmp_path), input_path])
    assert result == 7


def test_expression_refactor_noop(tmp_path: Path) -> None:
    """Functional test."""
    input_path1 = get_data("testdata/dockstore-tool-md5sum.cwl")
    input_path2 = get_data("testdata/echo-tool-packed.cwl")
    result = expression_refactor([str(tmp_path), input_path1, input_path2])
    assert result == 0


@pytest.fixture(scope="session")
def cwl_v1_0_dir(
    tmp_path_factory: TempPathFactory,
) -> Generator[str, None, None]:
    """Download the CWL 1.0.2 specs and return a path to the directory."""
    tmp_path = tmp_path_factory.mktemp("cwl_v1_0_dir")
    with cast(
        "HTTPResponse",
        requests.get(
            "https://github.com/common-workflow-language/common-workflow-language/archive/v1.0.2.tar.gz",
            stream=True,
        ).raw,
    ) as specfileobj:
        tf = tarfile.open(fileobj=specfileobj)
        if sys.version_info > (3, 12):
            tf.extractall(path=tmp_path, filter="data")
        else:
            tf.extractall(path=tmp_path)
    yield str(tmp_path / "common-workflow-language-1.0.2")
    shutil.rmtree(os.path.join(tmp_path))
