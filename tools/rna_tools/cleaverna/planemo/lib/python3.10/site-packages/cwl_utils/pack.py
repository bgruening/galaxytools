# SPDX-License-Identifier: Apache-2.0
#  Copyright (c) 2021 Michael R. Crusoe
#  Copyright (c) 2020 Seven Bridges
#  See https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/LICENSE
"""
CWL document packing functions.

The link resolution is as follows:

We always have two components: the base and the link
If the link is a url or absolute path it is what is used to fetch the data.
If the link is a relative path it is combined with the base and that is what is
used to fetch data

From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py
"""

import logging
import sys
import urllib.parse
import urllib.request
from collections.abc import ItemsView
from typing import TYPE_CHECKING, Any, Optional, Union, cast

from packaging import version

from cwl_utils import schemadef, utils

if TYPE_CHECKING:
    from _collections_abc import dict_items

logger = logging.getLogger(__name__)


def get_inner_dict(
    cwl: dict[str, Any], path: list[dict[str, Any]]
) -> Optional[dict[str, Any]]:
    if len(path) == 0:
        return cwl

    if isinstance(cwl, dict):
        _v = cwl.get(path[0]["key"])
        if _v is not None:
            return get_inner_dict(_v, path[1:])

    elif isinstance(cwl, list):  # Going to assume this is a map expressed as list
        for _v in cwl:
            if isinstance(_v, dict):
                if _v.get(path[0]["key_field"]) == path[0]["key"]:
                    return get_inner_dict(_v, path[1:])

    return None


def pack_process(
    cwl: dict[str, Any],
    base_url: urllib.parse.ParseResult,
    cwl_version: str,
    parent_user_defined_types: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    cwl = listify_everything(cwl)
    cwl = normalize_sources(cwl)
    cwl, user_defined_types = load_schemadefs(cwl, base_url, parent_user_defined_types)
    cwl = resolve_schemadefs(cwl, base_url, user_defined_types)
    cwl = resolve_imports(cwl, base_url)
    cwl = resolve_steps(
        cwl,
        base_url,
        cwl.get("cwlVersion", cwl_version),
        user_defined_types,
    )
    cwl = add_missing_requirements(cwl)
    return cwl


def listify_everything(cwl: dict[str, Any]) -> dict[str, Any]:
    """
    Convert many CWL construct from their map to the list version.

    See https://www.commonwl.org/v1.1/Workflow.html#map
    """
    for port in ["inputs", "outputs"]:
        cwl[port] = utils.normalize_to_list(
            cwl.get(port, []), key_field="id", value_field="type"
        )

    cwl["requirements"] = utils.normalize_to_list(
        cwl.get("requirements", []), key_field="class", value_field=None
    )

    if cwl.get("class") != "Workflow":
        return cwl

    cwl["steps"] = utils.normalize_to_list(
        cwl.get("steps", []), key_field="id", value_field=None
    )

    for _, v in enumerate(cwl["steps"]):
        if isinstance(v, dict):
            v["in"] = utils.normalize_to_list(
                v.get("in", []), key_field="id", value_field="source"
            )

    return cwl


def normalize_sources(cwl: dict[str, Any]) -> dict[str, Any]:
    """Normalize the steps and output of a CWL Workflow."""
    if cwl.get("class") != "Workflow":
        return cwl

    for _step in cwl.get("steps", {}):
        if not isinstance(_step, dict):
            continue

        _inputs = _step.get("in", {})
        for k, _input in enumerate(_inputs):
            if isinstance(_input, str):
                _inputs[k] = _normalize(_input)
            elif isinstance(_input, dict):
                _src = _input.get("source")
                if isinstance(_src, str):
                    _input["source"] = _normalize(_input["source"])

    _outputs = cwl.get("outputs", {})
    for k, _output in enumerate(_outputs):
        if isinstance(_output, str):
            _outputs[k] = _normalize(_output)
        elif isinstance(_output, dict):
            _src = _output.get("outputSource")
            if isinstance(_src, str):
                _output["outputSource"] = _normalize(_output["outputSource"])

    return cwl


def _normalize(s: str) -> str:
    if s.startswith("#"):
        return s[1:]
    else:
        return s


def load_schemadefs(
    cwl: dict[str, Any],
    base_url: urllib.parse.ParseResult,
    parent_user_defined_types: Optional[dict[str, Any]] = None,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Internalize any SchemaDefRequirement, and remove it."""
    user_defined_types = schemadef.build_user_defined_type_dict(cwl, base_url)
    if parent_user_defined_types is not None:
        user_defined_types.update(parent_user_defined_types)

    cwl["requirements"] = [
        req
        for req in cwl.get("requirements", [])
        if req.get("class") != "SchemaDefRequirement"
    ]

    return cwl, user_defined_types


def resolve_schemadefs(
    cwl: dict[str, Any],
    base_url: urllib.parse.ParseResult,
    user_defined_types: dict[str, Any],
) -> dict[str, Any]:
    cwl = schemadef.inline_types(cwl, "inputs", base_url, user_defined_types)
    cwl = schemadef.inline_types(cwl, "outputs", base_url, user_defined_types)
    return cwl


def resolve_imports(cwl: Any, base_url: urllib.parse.ParseResult) -> Any:
    if isinstance(cwl, dict):
        itr: Union["dict_items[Any, Any]", ItemsView[Any, Any]] = cwl.items()
    elif isinstance(cwl, list):
        itr = cast(ItemsView[Any, Any], [(n, v) for n, v in enumerate(cwl)])
    else:
        return cwl

    for k, v in itr:
        if isinstance(v, dict):
            if len(v) == 1:
                _k = list(v.keys())[0]
                if _k in ["$import", "$include"]:
                    cwl[k], this_base_url = utils.load_linked_file(
                        base_url, v[_k], is_import=_k == "$import"
                    )

        cwl[k] = resolve_imports(cwl[k], base_url)

    return cwl


def resolve_steps(
    cwl: dict[str, Any],
    base_url: urllib.parse.ParseResult,
    cwl_version: str,
    parent_user_defined_types: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    """Load and pack all "run" sections of the workflow steps."""
    if isinstance(cwl, str):
        raise RuntimeError(f"{base_url.geturl()}: Expecting a process, found a string")

    if not isinstance(cwl, dict):
        return cwl

    if cwl.get("class") != "Workflow":
        return cwl

    for _, v in enumerate(cwl["steps"]):
        if isinstance(v, dict):
            sys.stderr.write(
                f"\n--\nRecursing into step {base_url.geturl()}:{v['id']}\n"
            )

            _run = v.get("run")
            if isinstance(_run, str):
                v["run"], new_base_url = utils.load_linked_file(
                    base_url, _run, is_import=True
                )
                v["run"] = pack_process(
                    v["run"],
                    new_base_url,
                    cwl.get("cwlVersion", cwl_version),
                )
            else:
                v["run"] = pack_process(
                    v["run"],
                    base_url,
                    cwl.get("cwlVersion", cwl_version),
                    parent_user_defined_types,
                )
            if "cwlVersion" in v["run"]:
                parent_version = version.parse(
                    cwl.get("cwlVersion", cwl_version).strip("v")
                )
                this_version = version.parse(v["run"]["cwlVersion"].strip("v"))
                if this_version > parent_version:
                    cwl["cwlVersion"] = v["run"]["cwlVersion"]
                    # not really enough, but hope for the best

    return cwl


def add_missing_requirements(cwl: dict[str, Any]) -> dict[str, Any]:
    """Due to packing, we may need to add a "SubworkflowFeatureRequirement"."""
    requirements = cwl.get("requirements", [])
    present = {req["class"] for req in requirements}

    def _add_req(_req_name: str) -> None:
        nonlocal requirements
        if _req_name not in present:
            requirements += [{"class": _req_name}]

    if cwl.get("class") == "Workflow":
        sub_workflow = False
        for step in cwl["steps"]:
            if step["run"]["class"] == "Workflow":
                sub_workflow = True
                break
        if sub_workflow:
            _add_req("SubworkflowFeatureRequirement")
    return cwl


def pack(cwl_path: str) -> dict[str, Any]:
    """Pack a CWL document at the given path."""
    sys.stderr.write(f"Packing {cwl_path}\n")
    file_path_url = urllib.parse.urlparse(cwl_path)

    cwl, full_url = cast(
        tuple[dict[str, Any], urllib.parse.ParseResult],
        utils.load_linked_file(base_url=file_path_url, link="", is_import=True),
    )
    if "$graph" in cwl:
        # assume already packed
        return cwl
    cwl = pack_process(cwl, full_url, cwl["cwlVersion"])

    return cwl
