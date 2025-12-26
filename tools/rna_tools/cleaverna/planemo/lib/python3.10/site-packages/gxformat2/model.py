"""Abstractions for dealing with Format2 data."""
import logging
import os
from typing import (
    Any,
    cast,
    Optional,
    Union,
)

from typing_extensions import Literal

log = logging.getLogger(__name__)

DictOrList = Union[dict, list]
ConnectDict = dict


NativeGalaxyStepType = Literal[
    "subworkflow",
    "data_input",
    "data_collection_input",
    "tool",
    "pause",
    "parameter_input",
]
GxFormat2StepTypeAlias = Literal[
    "input",
    "input_collection",
    "parameter",
]
StepTypes = Union[NativeGalaxyStepType, GxFormat2StepTypeAlias]


STEP_TYPES = [
    "subworkflow",
    "data_input",
    "data_collection_input",
    "tool",
    "pause",
    "parameter_input",
]
STEP_TYPE_ALIASES: dict[GxFormat2StepTypeAlias, NativeGalaxyStepType] = {
    'input': 'data_input',
    'input_collection': 'data_collection_input',
    'parameter': 'parameter_input',
}


def get_native_step_type(gxformat2_step_dict: dict) -> NativeGalaxyStepType:
    """Infer native galaxy step type from the gxformat2 step as a dict."""
    specifies_subworkflow_run = bool(gxformat2_step_dict.get("run"))
    step_type_default = "tool" if not specifies_subworkflow_run else "subworkflow"
    raw_step_type = gxformat2_step_dict.get("type", step_type_default)
    if raw_step_type not in STEP_TYPES and raw_step_type not in STEP_TYPE_ALIASES:
        raise Exception(f"Unknown step type encountered {raw_step_type}")
    step_type:  NativeGalaxyStepType
    if raw_step_type in STEP_TYPE_ALIASES:
        step_type = STEP_TYPE_ALIASES[cast(GxFormat2StepTypeAlias, raw_step_type)]
    else:
        step_type = cast(NativeGalaxyStepType, raw_step_type)
    return step_type


# source: step#output and $link: step#output instead of outputSource: step/output and $link: step/output
SUPPORT_LEGACY_CONNECTIONS = os.environ.get("GXFORMAT2_SUPPORT_LEGACY_CONNECTIONS") == "1"


def pop_connect_from_step_dict(step: dict) -> ConnectDict:
    """Merge 'in' and 'connect' keys into a unified connection dict separated from state.

    Meant to be used an initial processing step in reasoning about connections defined by the
    format2 step description.
    """
    if "connect" not in step:
        step["connect"] = {}

    connect = step["connect"]
    del step["connect"]

    # handle CWL-style in dict connections.
    if "in" in step:
        step_in = step["in"]
        assert isinstance(step_in, dict)
        connection_keys = set()
        for key, value in step_in.items():
            # TODO: this can be a list right?
            if isinstance(value, dict) and 'source' in value:
                value = value["source"]
            elif isinstance(value, dict) and 'default' in value:
                continue
            elif isinstance(value, dict):
                raise KeyError(f'step input must define either source or default {value}')
            connect[key] = [value]
            connection_keys.add(key)

        for key in connection_keys:
            del step_in[key]

        if len(step_in) == 0:
            del step['in']

    return connect


def setup_connected_values(value, key: str = "", append_to: Optional[dict[str, list]] = None) -> Any:
    """Replace links with connected value."""

    def append_link(key: str, value: dict):
        if append_to is None:
            return

        if key not in append_to:
            append_to[key] = []

        assert "$link" in value
        link_value = value["$link"]
        append_to[key].append(clean_connection(link_value))

    def recurse(sub_value, sub_key) -> Any:
        return setup_connected_values(sub_value, sub_key, append_to=append_to)

    if _is_link(value):
        append_link(key, value)
        # Filled in by the connection, so to force late
        # validation of the field just mark as ConnectedValue,
        # which should be further validated by Galaxy
        return _connected_value()
    if isinstance(value, dict):
        new_dict_values: dict[str, Any] = {}
        for dict_k, dict_v in value.items():
            new_key = _join_prefix(key, dict_k)
            new_dict_values[dict_k] = recurse(dict_v, new_key)
        return new_dict_values
    elif isinstance(value, list):
        new_list_values: list[Any] = []
        for i, list_v in enumerate(value):
            # If we are a repeat we need to modify the key
            # but not if values are actually $links.
            if _is_link(list_v):
                assert isinstance(list_v, dict)
                append_link(key, list_v)
                new_list_values.append(None)
            else:
                new_key = "%s_%d" % (key, i)
                new_list_values.append(recurse(list_v, new_key))
        return new_list_values
    else:
        return value


def clean_connection(value: str) -> str:
    """Convert legacy style connection targets with modern CWL-style ones."""
    if value and "#" in value and SUPPORT_LEGACY_CONNECTIONS:
        # Hope these are just used by Galaxy testing workflows and such, and not in production workflows.
        log.warn(f"Legacy workflow syntax for connections [{value}] will not be supported in the future")
        value = value.replace("#", "/", 1)

    return value


def _connected_value():
    return {"__class__": "ConnectedValue"}


def _is_link(value: Any) -> bool:
    return isinstance(value, dict) and "$link" in value


def _join_prefix(prefix: Optional[str], key: str):
    if prefix:
        new_key = f"{prefix}|{key}"
    else:
        new_key = key
    return new_key


def convert_dict_to_id_list_if_needed(
    dict_or_list: DictOrList,
    add_label: bool = False,
    mutate: bool = False,
) -> list:
    """Convert a list or dict to a list with keys embedded.

    If `add_label` is True, embed dict keys as 'label' attribute
    else 'id'.
    """
    if isinstance(dict_or_list, dict):
        rval = []
        for key, value in dict_or_list.items():
            if not isinstance(value, dict):
                value = {"type": value}
            if not mutate:
                value = value.copy()
            if add_label:
                if value.get("label") is None:
                    value["label"] = key
            else:
                value["id"] = key
            rval.append(value)
    else:
        rval = cast(list, dict_or_list)
    return rval


def with_step_ids(steps: list, inputs_offset: int = 0):
    """Walk over a list of steps and ensure the steps have a numeric id if otherwise missing."""
    assert isinstance(steps, list)
    new_steps = []
    for i, step in enumerate(steps):
        if "id" not in step:
            step = step.copy()
            step["id"] = i + inputs_offset
        assert step["id"] is not None
        new_steps.append(step)
    return new_steps


def ensure_step_position(step: dict, order_index: int):
    """Ensure step contains a position definition.

    Modifies the input step dictionary.
    """
    if "position" not in step:
        step["position"] = {
            "left": 10 * order_index,
            "top": 10 * order_index
        }


def prune_position(step):
    """Keep only ``left`` and ``top`` keys in step position."""
    return {k: v for k, v in step.get('position', {}).items() if k in ('left', 'top')}


def native_input_to_format2_type(step: dict, tool_state: dict) -> Union[str, list[str]]:
    """Return a Format2 input type ('type') from a native input step dictionary."""
    module_type = step.get("type")
    if module_type == 'data_collection_input':
        format2_type = 'collection'
    elif module_type == 'data_input':
        format2_type = 'data'
    elif module_type == "parameter_input":
        native_type = cast(str, tool_state.get("parameter_type"))
        format2_type = native_type
        if native_type == "integer":
            format2_type = "int"
        elif native_type == "text":
            format2_type = "string"
        if tool_state.get("multiple", False):
            return [format2_type]
    return format2_type


def inputs_as_normalized_steps(workflow_dict):
    """Return workflow inputs to a steps in array.

    Normalize Format2 inputs. `workflow_dict` is a Format 2 representation of
    a workflow. This method does not modify `workflow_dict`.
    """
    if "inputs" not in workflow_dict:
        return []

    inputs = workflow_dict.get("inputs", [])
    new_steps = []
    inputs = convert_dict_to_id_list_if_needed(inputs)
    for input_def_raw in with_step_ids(inputs):
        input_def = input_def_raw.copy()

        if "label" in input_def and "id" in input_def:
            raise Exception("label and id are aliases for inputs, may only define one")
        if "label" not in input_def and "id" not in input_def:
            raise Exception("Input must define a label.")

        raw_label = input_def.pop("label", None)
        raw_id = input_def.pop("id", None)
        label = raw_label or raw_id

        if label is None:
            raise Exception("Input label must not be empty.")

        step_type = input_def.pop("type", "data")
        if step_type == "File":
            step_type = "data"
        elif step_type == "integer":
            step_type = "int"
        elif step_type == "text":
            step_type = "string"

        step_def = input_def
        step_def.update({
            "type": step_type,
            "id": label,
        })
        new_steps.append(step_def)

    return new_steps


def inputs_as_native_steps(workflow_dict: dict):
    """Return workflow inputs to a steps in array - like in native Galaxy.

    Convert Format2 types into native ones. `workflow_dict` is a Format 2
    representation of a workflow. This method does not modify `workflow_dict`.
    """
    if "inputs" not in workflow_dict:
        return []

    inputs = workflow_dict.get("inputs", [])
    new_steps = []
    inputs = convert_dict_to_id_list_if_needed(inputs)
    for input_def_raw in inputs:
        input_def = input_def_raw.copy()

        if "label" in input_def and "id" in input_def:
            raise Exception("label and id are aliases for inputs, may only define one")
        if "label" not in input_def and "id" not in input_def:
            raise Exception("Input must define a label.")

        raw_label = input_def.pop("label", None)
        raw_id = input_def.pop("id", None)
        label = raw_label or raw_id

        if label is None:
            raise Exception("Input label must not be empty.")

        input_type = input_def.pop("type", "data")
        if isinstance(input_type, list):
            if len(input_type) != 1:
                raise Exception("Only simple arrays of workflow inputs are currently supported")
            input_type = input_type[0]
            if input_type in ["File", "data", "data_input"]:
                raise Exception(f"Array of {input_type} is not supported")
            input_def["tool_state"] = {"multiple": True}
        if input_type in ["File", "data", "data_input"]:
            step_type = "data_input"
        elif input_type in ["collection", "data_collection", "data_collection_input"]:
            step_type = "data_collection_input"
        elif input_type in ["text", "string", "integer", "int", "float", "color", "boolean"]:
            step_type = "parameter_input"
            format2_type = input_type
            if format2_type == "int":
                native_type = "integer"
            elif format2_type == "string":
                native_type = "text"
            else:
                native_type = format2_type
            input_def["parameter_type"] = native_type
        else:
            raise Exception(f"Unknown input type [{input_type}] encountered.")

        step_def = input_def
        step_def.update({
            "type": step_type,
            "label": label,
        })
        default = step_def.get("default")
        if isinstance(default, dict) and default.get('class') == 'File':
            # First 'default' is input name, hardcoded to default, second 'default'
            # is the actual default for the input name
            step_def['in'] = {'default': {'default': step_def.pop('default')}}
        new_steps.append(step_def)

    return new_steps


def outputs_as_list(as_python: dict) -> list:
    """Extract outputs from Format2 rep as list."""
    outputs = as_python.get("outputs", [])
    outputs = convert_dict_to_id_list_if_needed(outputs)
    return outputs


def steps_as_list(format2_workflow: dict, add_ids: bool = False, inputs_offset: int = 0, mutate: bool = False) -> list[dict[str, Any]]:
    """Return steps as a list, converting ID map to list representation if needed.

    This method does mutate the supplied steps, try to make progress toward not doing this.

    Add keys as labels instead of IDs. Why am I doing this?
    """
    if "steps" not in format2_workflow:
        raise Exception(f"No 'steps' key in dict, keys are {format2_workflow.keys()}")
    steps = format2_workflow["steps"]
    steps = convert_dict_to_id_list_if_needed(steps, add_label=True, mutate=mutate)
    if add_ids:
        if mutate:
            append_step_id_to_step_list_elements(steps, inputs_offset=inputs_offset)
        else:
            steps = with_step_ids(steps, inputs_offset=inputs_offset)
    return steps


def append_step_id_to_step_list_elements(steps: list[dict[str, Any]], inputs_offset: int = 0) -> None:
    """Ensure a list of steps each contains an 'id' element."""
    assert isinstance(steps, list)
    for i, step in enumerate(steps):
        if "id" not in step:
            step["id"] = i + inputs_offset
        assert step["id"] is not None
