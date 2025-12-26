# SPDX-License-Identifier: Apache-2.0
"""Test that scoping of identifiers in Workflow.steps[].run is correct."""


from cwl_utils.parser import Workflow, load_document_by_uri

from .util import get_path


def test_workflow_step_process_scope_v1_0() -> None:
    """CWL v1.0 IDs under Workflow.steps[].run should not be scoped in the "run" scope."""
    uri = get_path("testdata/workflow_input_format_expr.cwl").as_uri()
    cwl_obj: Workflow = load_document_by_uri(uri)
    assert cwl_obj.steps[0].run.inputs[0].id.endswith("#format_extract/target")


def test_workflow_step_process_scope_v1_1() -> None:
    """CWL v1.1 IDs under Workflow.steps[].run should be scoped in the "run" scope."""
    uri = get_path("testdata/workflow_input_format_expr_v1_1.cwl").as_uri()
    cwl_obj: Workflow = load_document_by_uri(uri)
    assert cwl_obj.steps[0].run.inputs[0].id.endswith("#format_extract/run/target")


def test_workflow_step_process_scope_v1_2() -> None:
    """CWL v1.2 IDs under Workflow.steps[].run should be scoped in the "run" scope."""
    uri = get_path("testdata/workflow_input_format_expr_v1_2.cwl").as_uri()
    cwl_obj: Workflow = load_document_by_uri(uri)
    assert cwl_obj.steps[0].run.inputs[0].id.endswith("#format_extract/run/target")
