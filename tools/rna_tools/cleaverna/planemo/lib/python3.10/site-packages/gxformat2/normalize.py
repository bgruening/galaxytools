"""Abstractions for uniform across formats."""
import copy
from typing import Union

from gxformat2._scripts import ensure_format2
from gxformat2.converter import (
    steps_as_list,
)
from gxformat2.model import (
    inputs_as_normalized_steps,
    outputs_as_list,
)
from gxformat2.yaml import ordered_load

NON_INPUT_TYPES = ["tool", "subworkflow", "pause"]


class Inputs:
    """An abstraction around a Galaxy workflow's inputs."""

    def __init__(self, workflow_dict):
        """Parse normalized inputs from dictified workflow representation."""
        self._inputs = inputs_normalized(workflow_dict=workflow_dict)

    def is_an_input(self, target_label):
        """Return true if step label/id is a input.

        This does assume all labels that are numeric strings refer to the Galaxy 'id'
        (i.e. order_index). We should make sure Galaxy doesn't allow purely numeric
        labels so this doesn't get confusing and we can ensure correctness.
        """
        if isinstance(target_label, str):
            for input_def in self._inputs:
                label = input_def.get("label") or input_def.get("id")
                if target_label == label:
                    return True

            for input_def in self._inputs:
                label = input_def.get("label") or input_def.get("id")
                if label.isdigit() and target_label == int(label):
                    return True

        if isinstance(target_label, int):
            for index, input_def in enumerate(self._inputs):
                if target_label == index:
                    return True
                input_id = input_def.get("id")
                if isinstance(input_id, str):
                    if input_id.isdigit() and int(input_id) == target_label:
                        return True
                else:
                    if target_label == input_id:
                        return True

        return False

    @property
    def count(self):
        """Return the number of declared inputs for the workflow."""
        return len(self._inputs)


class NormalizedWorkflow:
    """Present a view of a Format2 workflow that has been normalized.

    In a normalized view:
    - Steps are a list and ids/labels have been populated
    - references to anonymous outputs have been relaced with the correct anonymous label

    TODO:
    - input: <type> has been converted to input: {type: <type>}
    - step inputs have been replaced with input declarations
    - normalize out declarations to dicts?
    """

    def __init__(self, input_workflow: dict):
        """Create a NormalizedWorkflow from supplied Format2 workflow."""
        self.input_workflow = input_workflow
        normalized_workflow = copy.deepcopy(input_workflow)
        _replace_anonymous_output_references(normalized_workflow)
        _ensure_implicit_step_outs(normalized_workflow)
        inputs = normalized_workflow.get("inputs", None)
        if inputs:
            normalized_workflow["inputs"] = inputs_as_normalized_steps(normalized_workflow)
        self.normalized_workflow_dict = normalized_workflow


def steps_normalized(workflow_dict=None, workflow_path=None):
    """Walk over a normalized step rep. across workflow formats."""
    workflow_dict = _ensure_format2(workflow_dict=workflow_dict, workflow_path=workflow_path)
    steps = steps_as_list(workflow_dict)
    return inputs_as_normalized_steps(workflow_dict) + steps


def inputs_normalized(**kwd):
    """Call steps_normalized and retain just the input steps normalized."""
    steps = steps_normalized(**kwd)
    input_steps = []
    for step in steps:
        step_type = step.get("type") or 'tool'
        if step_type in NON_INPUT_TYPES:
            continue

        input_steps.append(step)

    return input_steps


def outputs_normalized(**kwd):
    """Ensure Format2 and return outputs.

    Probably should go farther and normalize source -> outputSource,
    but doesn't yet do this.
    """
    workflow_dict = _ensure_format2(**kwd)
    return outputs_as_list(workflow_dict)


def walk_id_list_or_dict(dict_or_list: Union[dict, list]):
    """Walk over idmap regardless of list or dict representation."""
    if isinstance(dict_or_list, list):
        for item in dict_or_list:
            yield item["id"], item
    else:
        for item in dict_or_list.items():
            yield item


def _replace_anonymous_output_references(workflow_dict: dict):
    """Replace anonymous output referneces in a Format 2 representation.

    Don't do this when converting to/from other Galaxy formats because it
    effectively adds output labels, but it can be useful for checking logic
    and exporting to other formats that would benefit from or require an
    output label (e.g. abstract CWL).
    """
    runs_by_label = {}
    for step in steps_as_list(workflow_dict, add_ids=True, inputs_offset=len(workflow_dict["inputs"]), mutate=False):
        label = step.get("label")
        if label is None:
            label = step.get("id")
        assert label is not None, step
        label = str(label)
        if "run" in step:
            runs_by_label[label] = step["run"]

    for output_name, output in walk_id_list_or_dict(workflow_dict.get("outputs", {})):
        if "outputSource" in output:
            output_source = output["outputSource"]
            if "/" in output_source:
                step, output_name = output_source.split("/", 1)
                if ":" in output_name:
                    subworkflow_label, subworkflow_output = output_name.split(":", 1)
                    assert subworkflow_label in runs_by_label, f"{subworkflow_label} not in {runs_by_label.keys()}"
                    run = runs_by_label[subworkflow_label]
                    subworkflow_outputs = run["outputs"]
                    assert isinstance(subworkflow_outputs, dict)
                    for subworkflow_output_name, output_def in subworkflow_outputs.items():
                        if output_def["outputSource"] == f"{subworkflow_label}/{subworkflow_output}":
                            output["outputSource"] = f"{step}/{subworkflow_output_name}"


def _ensure_implicit_step_outs(workflow_dict: dict):
    """Ensure implicit 'out' dicts allowed by format2 are filled in for CWL."""
    outputs_by_label = {}

    def register_step_output(step_label, output_name):
        step_label = str(step_label)
        if step_label not in outputs_by_label:
            outputs_by_label[step_label] = set()
        outputs_by_label[step_label].add(output_name)

    def register_output_source(output_source):
        if "/" in output_source:
            step, output_name = output_source.split("/", 1)
            register_step_output(step, output_name)

    for output_name, output in walk_id_list_or_dict(workflow_dict.get("outputs", {})):
        if "outputSource" in output:
            output_source = output["outputSource"]
            if "/" in output_source:
                step, output_name = output_source.split("/", 1)
                register_step_output(step, output_name)

    for step in steps_as_list(workflow_dict, mutate=False):
        step_in = step.get("in", {})
        for step_in_def in step_in.values():
            register_output_source(step_in_def)

    for step in steps_as_list(workflow_dict, add_ids=True, inputs_offset=len(workflow_dict["inputs"]), mutate=True):
        label = step.get("label")
        if label is None:
            label = step.get("id")
        assert label is not None, step
        label = str(label)
        if "out" not in step:
            step["out"] = []
        for out in outputs_by_label.get(label, []):
            step_out = step["out"]
            if isinstance(step_out, list):
                if out not in step_out:
                    step_out.append(out)
            else:
                step_out[out] = {}


def _ensure_format2(workflow_dict=None, workflow_path=None):
    if workflow_path is not None:
        assert workflow_dict is None
        with open(workflow_path) as f:
            workflow_dict = ordered_load(f)

    workflow_dict = ensure_format2(workflow_dict)
    return workflow_dict
