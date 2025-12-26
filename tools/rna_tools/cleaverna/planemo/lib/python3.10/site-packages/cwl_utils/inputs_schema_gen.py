#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2024 Hirotaka Suetake
# Copyright 2024 Alexis Lucattini

"""Generate JSON Schema from CWL inputs object."""
import argparse
import json
import logging
import sys
from copy import deepcopy
from pathlib import Path
from typing import Any, Optional, Union
from urllib.parse import urlparse

import requests

# Get typeguard from extensions if we're running in python3.8
if sys.version_info[:2] < (3, 10):
    from typing_extensions import TypeGuard  # Not in 3.8 typing module
else:
    from typing import TypeGuard

from cwl_utils.loghandler import _logger as _cwlutilslogger
from cwl_utils.parser import (
    CommandLineTool,
    Directory,
    File,
    InputArraySchema,
    InputArraySchemaTypes,
    InputEnumSchema,
    InputEnumSchemaTypes,
    InputRecordSchema,
    InputRecordSchemaTypes,
    Workflow,
    WorkflowInputParameter,
    load_document_by_uri,
)
from cwl_utils.utils import (
    get_value_from_uri,
    is_local_uri,
    is_uri,
    sanitise_schema_field,
    to_pascal_case,
)

_logger = logging.getLogger("cwl-inputs-schema-gen")  # pylint: disable=invalid-name
defaultStreamHandler = logging.StreamHandler()  # pylint: disable=invalid-name
_logger.addHandler(defaultStreamHandler)
_logger.setLevel(logging.INFO)
_cwlutilslogger.setLevel(100)

# Globals

# Maps CWL types to JSON Schema types
PRIMITIVE_TYPES_MAPPING = {
    "boolean": "boolean",
    "string": "string",
    "int": "integer",
    "float": "number",
    "long": "number",
    "double": "number",
    "null": "null",
}

JSON_TEMPLATE_PATH = (
    Path(__file__)
    .parent.joinpath("./templates/workflow_input_json_schema_template.json")
    .absolute()
    .resolve()
)

# Some type hinting
InputType = Union[
    InputArraySchema, InputEnumSchema, InputRecordSchema, str, File, Directory
]


# Don't need type checking at runtime


class JSONSchemaProperty:
    """Generate a JSON schema property from a CWL input parameter."""

    def __init__(
        self,
        name: str,
        type_: Union[InputType, list[InputType], str, Any],
        description: Optional[str] = "",
        required: Optional[bool] = False,
    ):
        """Initialise the JSONSchemaProperty object."""
        # Initialise values
        self.name: str = name
        self.type_: Union[InputType, list[InputType], str, Any] = type_
        self.description = description
        self.required = required
        self.type_dict = self.generate_type_dict()

    def generate_type_dict(self) -> dict[str, Any]:
        """Generate the type dict for a property from a CWL input parameter type."""
        # If the type is a list and contains null, then the property is not required
        if isinstance(self.type_, list) and "null" in self.type_:
            self.required = False
            self.type_ = list(filter(lambda type_item: type_item != "null", self.type_))

            # Check if we're down to one item, we can then squeeze
            if len(self.type_) == 1:
                self.type_ = self.type_[0]

        # type_ is still a list therefore we offer multiple input types for this parameter
        if isinstance(self.type_, list):
            # We use the oneOf keyword to specify multiple types
            type_dict = self.generate_type_dict_from_type_list(self.type_)
        # type_ is a single type
        else:
            type_dict = self.generate_type_dict_from_type(self.type_)

        # Add in the description to the type dict
        type_dict.update({"description": self.description})

        return type_dict

    def generate_type_dict_from_type(self, type_item: Any) -> dict[str, Any]:
        """
        Generate the type dict for a property from a CWL input parameter type.

        We call this function for each type in the type_ list
        In the case there are multiple types, each dict is added to the oneOf list
        """
        # Primitive types should have a 1-1 mapping
        # Between an CWL Input Parameter type and a JSON schema type
        if isinstance(type_item, str):
            if type_item in PRIMITIVE_TYPES_MAPPING.keys():
                return {"type": PRIMITIVE_TYPES_MAPPING[type_item]}
            elif type_item in ["stdin"]:
                return {"$ref": "#/definitions/File"}
            elif type_item in ["File", "Directory", "Any"]:
                return {"$ref": f"#/definitions/{type_item}"}
            # When item is a record schema type
            elif is_uri(type_item):
                return {
                    "$ref": f"#/definitions/{to_pascal_case(get_value_from_uri(type_item))}"
                }
            else:
                raise ValueError(f"Unknown type: {type_item}")
        elif isinstance(type_item, InputArraySchemaTypes):
            return {
                "type": "array",
                "items": self.generate_type_dict_from_type(type_item.items),
            }
        elif isinstance(type_item, InputEnumSchemaTypes):
            return {
                "type": "string",
                "enum": list(
                    map(
                        lambda symbol_iter: get_value_from_uri(symbol_iter),
                        type_item.symbols,
                    )
                ),
            }
        elif isinstance(type_item, InputRecordSchemaTypes):
            if type_item.fields is None:
                return {"type": "object"}
            if not isinstance(type_item.fields, list):
                _cwlutilslogger.error(
                    "Expected fields of InputRecordSchemaType to be a list"
                )
                raise TypeError
            return {
                "type": "object",
                "properties": {
                    get_value_from_uri(prop.name): self.generate_type_dict_from_type(
                        prop.type_
                    )
                    for prop in type_item.fields
                },
            }
        elif isinstance(type_item, dict):
            # Nested import
            # {'$import': '../relative/path/to/schema'}
            if "$import" in type_item.keys():
                # This path is a relative path to import
                return {
                    "$ref": f"#/definitions/{to_pascal_case(get_value_from_uri(type_item['$import']))}"
                }
            else:
                raise ValueError(f"Unknown type: {type_item}")
        elif isinstance(type_item, list):
            # Nested schema
            return {
                "oneOf": list(
                    map(
                        lambda type_iter: self.generate_type_dict_from_type(type_iter),
                        type_item,
                    )
                )
            }
        else:
            raise ValueError(f"Unknown type: {type_item}")

    def generate_type_dict_from_type_list(
        self, type_: list[InputType]
    ) -> dict[str, Any]:
        """Given a list of types, generate a JSON schema property dict wrapped in oneOf list."""
        return {
            "oneOf": list(
                map(
                    lambda type_item: self.generate_type_dict_from_type(type_item),
                    type_,
                )
            )
        }

    def to_dict(self) -> dict[str, Any]:
        """Return as a dictionary."""
        return {self.name: self.type_dict}


def get_is_required_from_input_parameter(
    input_parameter: WorkflowInputParameter,
) -> bool:
    """Given an input parameter, return if it is required."""
    if isinstance(input_parameter.type_, str) and input_parameter.type_.endswith("?"):
        return False
    if input_parameter.default is not None:
        return False
    if isinstance(input_parameter.type_, list) and "null" in input_parameter.type_:
        return False
    if isinstance(input_parameter.type_, InputRecordSchemaTypes):
        if input_parameter.type_ is not None:
            if (isinstance(input_parameter.type_.type_, str)) and (
                input_parameter.type_.type_.endswith("?")
            ):
                return False
    return True


def generate_json_schema_property_from_input_parameter(
    input_parameter: WorkflowInputParameter,
) -> JSONSchemaProperty:
    """
    Given an input parameter, generate a JSON schema property.

    :param input_parameter:
    :return:
    """
    # Get the input name and documentation for description
    input_name = get_value_from_uri(str(input_parameter.id))
    doc = input_parameter.doc
    required = get_is_required_from_input_parameter(input_parameter)

    return JSONSchemaProperty(
        name=input_name,
        type_=input_parameter.type_,
        description=doc if doc is not None else "",
        required=required,
    )


def generate_definition_from_schema(schema: InputRecordSchema) -> dict[str, Any]:
    """
    Given a schema, generate a JSON schema definition.

    :param schema:
    :return:
    """
    # Sanitise each field of the schema
    sanitised_fields = {}

    if schema.fields is None:
        return {}

    for field in schema.fields:
        sanitised_fields.update(
            {
                get_value_from_uri(field.name): sanitise_schema_field(
                    {"type": field.type_}
                )
            }
        )

    # Generate JSON properties
    property_list = []

    for prop_name, prop_obj in sanitised_fields.items():
        # Simplify type first by removing nulls
        required = True

        # If the property object is a string, then it's a reference to another schema
        if isinstance(prop_obj, str):
            raise TypeError("Property Object should be a dictionary")

        if isinstance(prop_obj.get("type", []), list):
            if "null" in prop_obj.get("type", []):
                required = False
            prop_obj["type"] = list(
                filter(lambda type_item: type_item != "null", prop_obj.get("type", []))
            )

            # Check if we're down to one item
            if len(prop_obj["type"]) == 1:
                prop_obj["type"] = prop_obj["type"][0]

        # Generate JSONSchema Property
        prop = JSONSchemaProperty(
            name=prop_name,
            type_=prop_obj.get("type"),
            description=prop_obj.get("doc", ""),
            required=required,
        )
        property_list.append(prop)

    return {
        to_pascal_case(get_value_from_uri(str(schema.name))): {
            "type": "object",
            "properties": {prop.name: prop.type_dict for prop in property_list},
            "required": [prop.name for prop in property_list if prop.required],
        }
    }


def cwl_to_jsonschema(cwl_obj: Union[Workflow, CommandLineTool]) -> Any:
    """
    cwl_obj: A CWL Object.

    Returns:
        A JSONSchema object.

    Example:
        cwl_obj = load_document_by_uri(<CWL_URL>)
        jsonschema = cwl_to_jsonschema(cwl_inputs)

    """
    # Initialise the schema from the workflow input json schema template
    with open(JSON_TEMPLATE_PATH) as template_h:
        input_json_schema = json.load(template_h)

    # Get the complex schema keys
    def is_complex_record_schema_key(idx_iter: str) -> TypeGuard[bool]:
        if cwl_obj.loadingOptions.idx is None:
            return False

        if cwl_obj.loadingOptions.idx.get(idx_iter) is None:
            return False

        if not isinstance(cwl_obj.loadingOptions.idx.get(idx_iter), tuple):
            return False

        # Get index as a tuple
        input_schema_type, _ = cwl_obj.loadingOptions.idx.get(idx_iter, (None, None))

        if isinstance(input_schema_type, InputRecordSchemaTypes):
            return True
        return False

    complex_schema_keys: list[str] = list(
        filter(
            lambda idx_iter: is_complex_record_schema_key(idx_iter),
            cwl_obj.loadingOptions.idx.keys(),
        )
    )

    # Complex schema values
    def get_complex_schema_values(idx_iter: str) -> InputRecordSchema:
        if not isinstance(cwl_obj.loadingOptions.idx.get(idx_iter), tuple):
            raise TypeError(f"Expected tuple from idx loading options key {idx_iter}")

        # Collect input record schema
        input_record_schema, _ = cwl_obj.loadingOptions.idx.get(idx_iter, (None, None))

        if not isinstance(input_record_schema, InputRecordSchemaTypes):
            raise TypeError(
                f"Expected InputRecordSchemaTypes from idx loading options key {idx_iter}"
            )

        return input_record_schema

    complex_schema_values: list[InputRecordSchema] = list(
        map(
            lambda idx_iter: get_complex_schema_values(idx_iter),
            complex_schema_keys,
        )
    )

    # Load in all $imports to be referred by complex input types
    workflow_schema_definitions_list = list(
        map(
            lambda complex_schema_values_iter: generate_definition_from_schema(
                complex_schema_values_iter
            ),
            complex_schema_values,
        )
    )

    if cwl_obj.requirements is not None:
        try:
            schema_def_requirement = next(
                filter(
                    lambda requirement_iter: requirement_iter.class_
                    == "SchemaDefRequirement",
                    cwl_obj.requirements,
                )
            )

            workflow_schema_definitions_list.extend(
                list(
                    map(
                        lambda schema_def_iter: generate_definition_from_schema(
                            schema_def_iter
                        ),
                        schema_def_requirement.types,
                    )
                )
            )

        except StopIteration:
            pass

    # Convert schema definitions to dict
    workflow_schema_definitions_dict = {}
    for schema_definition in workflow_schema_definitions_list:
        workflow_schema_definitions_dict.update(schema_definition)

    # Generate JSON Schema Properties
    properties = list(
        map(
            lambda workflow_parameter_input_obj: generate_json_schema_property_from_input_parameter(
                workflow_parameter_input_obj
            ),
            cwl_obj.inputs,
        )
    )

    # Generate JSON schema
    input_json_schema.update(
        {
            "type": "object",
            "properties": {
                prop.name: (
                    {"oneOf": [{"type": "null"}, prop.type_dict]}
                    if prop.required is False
                    else prop.type_dict
                )
                for prop in properties
            },
            "required": [prop.name for prop in properties if prop.required],
        }
    )

    # Update definitions from schema
    input_json_schema["definitions"].update(workflow_schema_definitions_dict)

    # Slim down the schema as required
    input_json_schema = slim_definitions(input_json_schema)

    # Add "additionalProperties": false to top of schema
    # input_json_schema["additionalProperties"] = False

    return input_json_schema


# Traverse the properties and return all definitions that are used
def _recursive_search(
    json_data: dict[str, Any],
    target_key: str,
) -> list[Any]:
    """Given a target key return all instances of a key in a json object."""
    result = []

    if isinstance(json_data, dict):
        for key, value in json_data.items():
            if key == target_key:
                result.append(value)
            else:
                result.extend(_recursive_search(value, target_key))
    elif isinstance(json_data, list):
        for item in json_data:
            result.extend(_recursive_search(item, target_key))

    return result


# Get all the property dependencies
def _get_all_ref_attributes(json_object: dict[str, Any]) -> list[Any]:
    """Given a json object, return all the reference attributes."""
    return _recursive_search(json_object, "$ref")


def get_property_dependencies(
    property_dict: dict[str, Any],
    input_json_schema: dict[str, Any],
    existing_property_dependencies: Optional[list[Any]] = None,
) -> list[str]:
    """Recursively collect all dependencies for a property."""
    # Initialise return list
    if existing_property_dependencies is None:
        existing_property_dependencies = []

    # All reference attributes
    for reference_attribute in _get_all_ref_attributes(property_dict):
        # Get the value from the reference attribute
        reference_value = get_value_from_uri(reference_attribute)
        # If the reference value is not in the existing property dependencies, add it
        if reference_value not in existing_property_dependencies:
            existing_property_dependencies.append(reference_value)
            # Get the property dependencies of the reference value
            existing_property_dependencies.extend(
                get_property_dependencies(
                    input_json_schema["definitions"][reference_value],
                    input_json_schema,
                    existing_property_dependencies,
                )
            )

    return existing_property_dependencies


def slim_definitions(input_json_schema: dict[str, Any]) -> dict[str, Any]:
    """
    Slim down the schema to only the definitions that are used by the properties.

    Traverse the properties and return all definitions that are used.
    Remove all other definitions.
    """
    # Copy schema
    input_json_schema = deepcopy(input_json_schema)

    # Get required definitions
    required_definitions = get_property_dependencies(
        input_json_schema.get("properties", {}), input_json_schema
    )

    for definition_key in list(input_json_schema["definitions"].keys()):
        if definition_key not in required_definitions:
            del input_json_schema["definitions"][definition_key]

    return input_json_schema


def arg_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(description="Generate JSON Schema from a CWL URI.")
    parser.add_argument("cwl_url", help="URL or Path to the CWL document")
    parser.add_argument(
        "-o",
        "--output",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Output file. Default is stdout.",
    )
    return parser


def parse_args(args: list[str]) -> argparse.Namespace:
    """Parse the command line arguments."""
    return arg_parser().parse_args(args)


def main() -> None:
    """Console entry point."""
    sys.exit(run(parse_args(sys.argv[1:])))


def get_cwl_url(url: str) -> str:
    """
    Conform to uri format.

    If no scheme, then assert is a local file path and exists
    if scheme is file, then assert is a local file path and exists
    If scheme is not file, then assert is a valid Web URL
    Return either the url or the local path as a uri.
    """
    if not is_uri(url):
        if not Path(url).exists():
            logging.error("The CWL URL is invalid.")
            raise FileNotFoundError
        return Path(url).as_uri()
    elif is_local_uri(url):
        if not Path(urlparse(url).path).exists():
            logging.error("The CWL URL is invalid.")
            raise FileNotFoundError
        return url
    else:
        # urlparse(url).scheme not in ['file']:
        response = requests.get(url, timeout=20)
        if response.status_code != 200:
            logging.error("The CWL URL is invalid.")
            raise FileNotFoundError
        return url


def run(args: argparse.Namespace) -> int:
    """Run the main program."""
    # Check the cwl_url is valid
    cwl_url = get_cwl_url(args.cwl_url)

    # Check the output file is writable
    if args.output.name != "<stdout>":
        if not Path(args.output.name).parent.is_dir():
            logging.error(
                "The output file is not writable, the output parent directory does not exist"
            )
            return 1

    _logger.info("Loading the CWL document")
    cwl_obj = load_document_by_uri(cwl_url)

    try:
        jsonschema = cwl_to_jsonschema(cwl_obj)
    except Exception as e:
        _logger.exception(
            "Failed to generate JSON Schema from CWL inputs object. Error: %s", e
        )
        return 1
    args.output.write(json.dumps(jsonschema, indent=2) + "\n")

    return 0


if __name__ == "__main__":
    main()
