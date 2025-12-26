# SPDX-License-Identifier: Apache-2.0
#  Copyright (c) 2023 Genomics plc
#  Copyright (c) 2021 Michael R. Crusoe
#  Copyright (c) 2020 Seven Bridges
#  See https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/LICENSE
"""
Valid forms of user defined types stored in external file.

A single dictionary (tests/types/singletype.yml)
A list of dictionaries (e.g. tests/types/recursive.yml)
Types can refer to other types in the file
Names can not clash across files (This seems arbitrary and we allow that for packing)
Only records and arrays can be defined (https://github.com/common-workflow-language/cwl-v1.2/pull/14)

From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py
"""


import sys
import urllib.parse
from copy import deepcopy
from typing import Any, cast

from cwl_utils import errors, types, utils


def build_user_defined_type_dict(
    cwl: dict[str, Any], base_url: urllib.parse.ParseResult
) -> dict[str, Any]:
    user_defined_types = {}
    # Check for `$import` directly under `requirements` so we can specially handle
    # the new base_url
    for index, entry in enumerate(cwl.get("requirements", [])):
        if isinstance(entry, dict) and "$import" in entry:
            requirement, new_base_url = utils.load_linked_file(
                base_url, entry["$import"], is_import=True
            )
            if requirement["class"] == "SchemaDefRequirement":
                cwl["requirements"][index] = requirement
                type_definition_list = requirement["types"]
                path_prefix = new_base_url.geturl()
                sys.stderr.write(
                    f"Parsing {len(type_definition_list)} types from {path_prefix}\n"
                )
                for v in type_definition_list:
                    k = v.get("name")
                    if k is None:
                        raise RuntimeError(f"In file {path_prefix} type missing name")
                    user_defined_types[f"{path_prefix}#{k}"] = v
    return _build_user_defined_type_dict(cwl, base_url, user_defined_types)


def _build_user_defined_type_dict(
    cwl: dict[str, Any],
    base_url: urllib.parse.ParseResult,
    user_defined_types: dict[str, Any],
) -> dict[str, Any]:
    schemadef: dict[str, str] = next(
        (
            req
            for req in cwl.get("requirements", [])
            if req.get("class") == "SchemaDefRequirement"
        ),
        {},
    )
    schema_list = cast(list[dict[str, Any]], schemadef.get("types", []))

    if not isinstance(schema_list, list):
        raise RuntimeError(
            f"In file {base_url.geturl()}: "
            f"Schemadef types have to be a list\n"
            f"Instead, got: {schema_list}"
        )

    for schema in schema_list:
        if not isinstance(schema, dict):
            raise RuntimeError(
                f"In file {base_url.geturl()}: "
                f"User type has to be a dict\n"
                f"Instead, got: {schema}"
            )

        if len(schema.keys()) == 1 and list(schema.keys())[0] == "$import":
            type_definition_list, this_url = utils.load_linked_file(
                base_url, schema["$import"], is_import=True
            )
            # This is always a list
            if isinstance(type_definition_list, dict):
                type_definition_list = [type_definition_list]
                # except when it isn't

            path_prefix = (
                this_url.geturl()
            )  # sbpack.lib.normalized_path(schema["$import"], base_url).geturl()
            sys.stderr.write(
                f"Parsing {len(type_definition_list)} types from {path_prefix}\n"
            )
            for v in type_definition_list:
                k = v.get("name")
                if k is None:
                    raise RuntimeError(f"In file {path_prefix} type missing name")
                user_defined_types[f"{path_prefix}#{k}"] = v

        else:
            path_prefix = base_url.geturl()
            user_defined_types[f"{path_prefix}#{schema.get('name')}"] = schema

    # sys.stderr.write(str(user_defined_types))
    # sys.stderr.write("\n")

    return user_defined_types


# port = "input" or "output"
def inline_types(
    cwl: dict[str, Any],
    port: str,
    base_url: urllib.parse.ParseResult,
    user_defined_types: dict[str, Any],
) -> dict[str, Any]:
    if (
        len(cwl[port]) == 1
        and isinstance(cwl[port][0], dict)
        and cwl[port][0]["id"] == "$import"
    ):
        defs, base_url = utils.load_linked_file(
            base_url, cwl[port][0]["type"], is_import=True
        )
    else:
        defs = cwl[port]

    cwl[port] = [_inline_type(v, base_url, user_defined_types) for v in defs]
    return cwl


_inline_type_name_uniq_id = 0
_inline_type_names: set[str] = set()


def _inline_type(
    v: Any, base_url: urllib.parse.ParseResult, user_defined_types: dict[str, Any]
) -> Any:
    global _inline_type_name_uniq_id

    _inline_type_name_uniq_id += 1

    if isinstance(v, str):
        # Handle syntactic sugar
        if v.endswith("[]"):
            return {
                "type": "array",
                "items": _inline_type(v[:-2], base_url, user_defined_types),
            }

        if v.endswith("?"):
            return ["null", _inline_type(v[:-1], base_url, user_defined_types)]

        if v in types.built_in_types:
            return v

        if "#" not in v:
            path_prefix = base_url
            path_suffix = v
        else:
            parts = v.split("#")
            path_prefix = utils.resolved_path(base_url, parts[0])
            path_suffix = parts[1]

        path = f"{path_prefix.geturl()}#{path_suffix}"

        if path not in user_defined_types:
            raise RuntimeError(
                f"Could not find type {path!r} in {user_defined_types!r}."
            )
        else:
            resolve_type = deepcopy(user_defined_types[path])
            # resolve_type.pop("name", None) # Should work, but cwltool complains
            if "name" in resolve_type:
                user_type_name = resolve_type["name"]
                if user_type_name in _inline_type_names:
                    resolve_type["name"] = (
                        f"{user_type_name}_{_inline_type_name_uniq_id}"
                    )
                else:
                    _inline_type_names.add(user_type_name)
            else:
                resolve_type["name"] = f"user_type_{_inline_type_name_uniq_id}"
            return _inline_type(resolve_type, path_prefix, user_defined_types)

    elif isinstance(v, list):
        return [_inline_type(_v, base_url, user_defined_types) for _v in v]

    elif isinstance(v, dict):
        if v.get("$import") is not None:
            imported_type, import_base_url = utils.load_linked_file(
                base_url, v["$import"], is_import=True
            )
            return _inline_type(imported_type, import_base_url, user_defined_types)

        _type = v.get("type")
        if _type is None:
            raise errors.MissingTypeName(
                f"In file {base_url.geturl()}, type {v.get('name')} is missing type name"
            )

        elif _type == "enum":
            return v

        elif _type == "array":
            if "items" not in v:
                raise errors.ArrayMissingItems(
                    f"In file {base_url.geturl()}, array type {_type.get('name')} is missing 'items'"
                )

            v["items"] = _inline_type(v["items"], base_url, user_defined_types)
            return v

        elif _type == "record":
            if "fields" not in v:
                raise errors.RecordMissingFields(
                    f"In file {base_url.geturl()}, record type {_type.get('name')} is missing 'fields'"
                )

            fields = utils.normalize_to_list(
                v["fields"], key_field="name", value_field="type"
            )
            v["fields"] = [
                _inline_type(_f, base_url, user_defined_types) for _f in fields
            ]
            return v

        elif _type in types.built_in_types:
            return v

        else:
            v["type"] = _inline_type(_type, base_url, user_defined_types)
            return v

    else:
        raise RuntimeError("Found a type sbpack can not understand")
