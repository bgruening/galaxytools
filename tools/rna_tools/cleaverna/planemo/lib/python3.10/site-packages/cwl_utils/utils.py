# SPDX-License-Identifier: Apache-2.0
"""Miscellaneous utility functions."""
import os
import pathlib
import subprocess  # nosec
import sys
import urllib.error
import urllib.parse
import urllib.request
from collections.abc import MutableMapping, MutableSequence
from copy import deepcopy
from io import StringIO
from typing import Any, Optional, Union
from urllib.parse import urlparse

from ruamel.yaml.main import YAML
from ruamel.yaml.parser import ParserError
from ruamel.yaml.scanner import ScannerError

from cwl_utils.errors import MissingKeyField
from cwl_utils.loghandler import _logger

# Type hinting
from cwl_utils.parser import InputRecordSchemaTypes

# Load as 1.2 files
from cwl_utils.parser.cwl_v1_2 import InputArraySchema as InputArraySchemaV1_2
from cwl_utils.parser.cwl_v1_2 import InputEnumSchema as InputEnumSchemaV1_2

fast_yaml = YAML(typ="safe")

_USERNS: Optional[bool] = None


def _is_github_symbolic_link(base_url: urllib.parse.ParseResult, contents: str) -> bool:
    """
    Test if link is a GitHub style symbolic link.

    Look for remote path with contents that is a single line with no new
    line with an extension.

    https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py
    """
    if base_url.scheme in ["file://", ""]:
        return False

    idx = contents.find("\n")
    if idx > -1:
        return False

    if "." not in contents:
        return False

    return True


def bytes2str_in_dicts(
    inp: Union[MutableMapping[str, Any], MutableSequence[Any], Any],
) -> Union[str, MutableSequence[Any], MutableMapping[str, Any]]:
    """
    Convert any present byte string to unicode string, inplace.

    input is a dict of nested dicts and lists
    """
    # if input is dict, recursively call for each value
    if isinstance(inp, MutableMapping):
        for k in inp:
            inp[k] = bytes2str_in_dicts(inp[k])
        return inp

    # if list, iterate through list and fn call
    # for all its elements
    if isinstance(inp, MutableSequence):
        for idx, value in enumerate(inp):
            inp[idx] = bytes2str_in_dicts(value)
            return inp

    # if value is bytes, return decoded string,
    elif isinstance(inp, bytes):
        return inp.decode("utf-8")

    # simply return elements itself
    return inp


def load_linked_file(
    base_url: urllib.parse.ParseResult, link: str, is_import: bool = False
) -> tuple[Any, urllib.parse.ParseResult]:
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""
    new_url = resolved_path(base_url, link)

    if new_url.scheme in ["file://", ""]:
        contents = pathlib.Path(new_url.path).read_text()
    else:
        try:
            contents = (
                urllib.request.urlopen(new_url.geturl()).read().decode("utf-8")  # nosec
            )
        except urllib.error.HTTPError as e:
            _logger.error("Could not find linked file: %s", new_url.geturl())
            raise SystemExit(e) from e

    if _is_github_symbolic_link(new_url, contents):
        # This is an exception for symbolic links on github
        sys.stderr.write(
            f"{new_url.geturl()}: found file-like string in contents.\n"
            f"Treating as github symbolic link to {contents}\n"
        )
        return load_linked_file(new_url, contents, is_import=is_import)

    if is_import:
        try:
            _node = fast_yaml.load(contents)
        except ParserError as e:
            e.context = f"\n===\nMalformed file: {new_url.geturl()}\n===\n" + e.context
            raise SystemExit(e) from e
        except ScannerError as e:
            e.problem = f"\n===\nMalformed file: {new_url.geturl()}\n===\n" + e.problem
            raise SystemExit(e) from e

    else:
        _node = contents

    return _node, new_url


def normalize_to_map(
    obj: Union[list[Any], dict[str, Any]], key_field: str
) -> dict[str, Any]:
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""
    if isinstance(obj, dict):
        return deepcopy(obj)
    elif isinstance(obj, list):
        map_obj = {}
        for v in obj:
            if not isinstance(v, dict):
                raise RuntimeError("Expecting a dict here")
            k = v.get(key_field)
            if k is None:
                raise MissingKeyField(key_field)
            v.pop(key_field, None)
            map_obj[k] = v
        return map_obj
    else:
        raise RuntimeError("Expecting a dictionary or a list here")


def normalize_to_list(
    obj: Union[list[Any], dict[str, Any]], key_field: str, value_field: Optional[str]
) -> list[Any]:
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""
    if isinstance(obj, list):
        return deepcopy(obj)
    elif isinstance(obj, dict):
        map_list = []
        for k, v in obj.items():
            if not isinstance(v, dict):
                if value_field is None:
                    raise RuntimeError(f"Expecting a dict here, got {v}")
                v = {value_field: v}
            v.update({key_field: k})
            map_list += [v]
        return map_list
    else:
        raise RuntimeError("Expecting a dictionary or a list here")


def resolved_path(
    base_url: urllib.parse.ParseResult, link: str
) -> urllib.parse.ParseResult:
    """
    Derive a resolved path.

    This function will
    1. Resolve the path, which means dot and double dot components are resolved
    2. Use the OS appropriate path resolution for local paths, and network
    appropriate resolution for network paths

    From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py

    :param base_url: "this document"
    :param link: "string in this document"
    :returns: new URL that allows us to retrieve the linked document
    """
    link_url = urllib.parse.urlparse(link)
    # The link will always Posix

    if link_url.scheme == "file://":
        # Absolute local path
        return link_url

    elif link_url.scheme == "":
        # Relative path, can be local or remote
        if base_url.scheme in ["file://", ""]:
            # Local relative path
            if link == "":
                return base_url
            else:
                return urllib.parse.urlparse(
                    urllib.parse.urljoin(base_url.geturl(), link_url.geturl())
                )
        else:
            # Remote relative path
            return urllib.parse.urlparse(
                urllib.parse.urljoin(base_url.geturl(), link_url.path)
            )
            # We need urljoin because we need to resolve relative links in a
            # platform independent manner

    # Absolute remote path
    return link_url


def singularity_supports_userns() -> bool:
    """Confirm if the version of Singularity install supports the --userns flag."""
    global _USERNS  # pylint: disable=global-statement
    if _USERNS is None:
        try:
            hello_image = os.path.join(os.path.dirname(__file__), "hello.simg")
            result = subprocess.Popen(  # nosec
                ["singularity", "exec", "--userns", hello_image, "true"],
                stderr=subprocess.PIPE,
                stdout=subprocess.DEVNULL,
                universal_newlines=True,
            ).communicate(timeout=60)[1]
            _USERNS = (
                "No valid /bin/sh" in result
                or "/bin/sh doesn't exist in container" in result
                or "executable file not found in" in result
            )
        except subprocess.TimeoutExpired:
            _USERNS = False
    return _USERNS


def yaml_dumps(obj: Any) -> str:
    """
    Shortcut.

    Don't use if you have a file descriptor (like sys.stdout) available.
    """
    yaml = YAML()
    stream = StringIO()
    yaml.dump(obj, stream)
    return stream.getvalue()


def to_pascal_case(name: str) -> str:
    """
    Convert a string to PascalCase.

    fastq-list-row to FastqListRow
    fastq_list_row to FastqListRow
    :param name:
    :return:
    """
    return "".join(
        map(lambda word: word.capitalize(), name.replace("_", "-").split("-"))
    )


def sanitise_schema_field(
    schema_field_item: Union[dict[str, Any], str],
) -> Union[dict[str, Any], str]:
    """
    Schemas need to be resolved before converted to JSON properties.

    Convert
      {
        'type': 'Directory?'
      }
    To
      {
        'type': ['null', 'Directory']
      }

    Convert
      {
        'type': 'string[]'
      }
    To
      InputArraySchema(
        type_=array,
        items=string
      )

    Convert
      {
        'type': 'File[]?'
      }
    To
      {
        'type': [
          'null', InputArraySchema(
            type_=array,
            items=File
          )
        ]
      }

    Convert
      {
        'type': 'Enum',
        'symbols': ['A', 'B', 'C']
      }

    To
      {
        'type': InputEnumSchema(
          type_=enum,
          symbols=['A', 'B', 'C']
        )
      }

    Convert
      {
        'type': 'array',
        'items': {
          '$import': '../../../schemas/fastq-list-row/1.0.0/fastq-list-row__1.0.0.yaml#fastq-list-row'
        }
      }
    To
      {
        'type': InputArraySchema(
          type_=array,
          items={
            '$import': '../../../schemas/fastq-list-row/1.0.0/fastq-list-row__1.0.0.yaml#fastq-list-row'
          }
        )
      }

    :param schema_field_item:
    :return:
    """
    # We might be just a string, in which case, just return
    # This happens in the case that type is a list of primitive types
    if isinstance(schema_field_item, str):
        return schema_field_item

    # Copy schema field
    schema_field_item = deepcopy(schema_field_item)
    required = True

    if isinstance(schema_field_item, InputRecordSchemaTypes):
        return schema_field_item

    if isinstance(schema_field_item.get("type"), list):
        if "null" in schema_field_item.get("type", []):
            required = False
        schema_field_item["type"] = list(
            filter(
                lambda type_item: type_item != "null", schema_field_item.get("type", [])
            )
        )
        if len(schema_field_item["type"]) == 1:
            schema_field_item["type"] = schema_field_item["type"][0]
        else:
            # Recursively get items
            schema_field_item["type"] = list(
                map(
                    lambda field_subtypes: sanitise_schema_field(field_subtypes),
                    schema_field_item.get("type", []),
                )
            )

    if isinstance(schema_field_item.get("type"), str):
        if schema_field_item.get("type", "").endswith("?"):
            required = False
            schema_field_item["type"] = schema_field_item.get("type", "").replace(
                "?", ""
            )

        if schema_field_item.get("type", "").endswith("[]"):
            # Strip list
            schema_field_item["type"] = schema_field_item.get("type", "").replace(
                "[]", ""
            )
            # Convert to array
            schema_field_item["type"] = InputArraySchemaV1_2(
                type_="array", items=schema_field_item.get("type", "")
            )

    if isinstance(schema_field_item.get("type"), dict):
        # Likely an enum
        if schema_field_item.get("type", {}).get("type", "") == "enum":
            schema_field_item["type"] = InputEnumSchemaV1_2(
                type_="enum",
                symbols=schema_field_item.get("type", {}).get("symbols", ""),
            )
        elif schema_field_item.get("type", {}).get("type", "") == "array":
            schema_field_item["type"] = InputArraySchemaV1_2(
                type_="array", items=schema_field_item.get("type", {}).get("items", "")
            )
        elif "$import" in schema_field_item.get("type", {}).keys():
            # Leave import as is
            pass
        else:
            raise ValueError(f"Unknown type: {schema_field_item.get('type')}")

    if not required:
        if isinstance(schema_field_item.get("type"), list):
            schema_field_item["type"] = ["null"] + schema_field_item.get("type", [])
        else:
            schema_field_item["type"] = ["null", schema_field_item.get("type", "")]

    return schema_field_item


def is_uri(uri: str) -> bool:
    """
    Given a URI return True if it is a URI.

    :param uri:
    :return:
    """
    if not urlparse(uri).scheme == "":
        return True
    else:
        return False


def is_local_uri(uri: str) -> bool:
    """Given a uri, first check if it is a uri, then check if it is a local uri."""
    if is_uri(uri) and urlparse(uri).scheme == "file":
        return True
    return False


def get_value_from_uri(uri: str) -> str:
    """
    Given a URI, return the value after #.

    file://path/to/imported/record#my_workflow_name/record_name
    Returns
    record_name
    :param uri:
    :return:
    """
    url_obj = urlparse(uri)
    return url_obj.fragment.rsplit("/")[-1]
