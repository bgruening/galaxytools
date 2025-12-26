# SPDX-License-Identifier: Apache-2.0
"""Tests for cwl-inputs-schema-gen."""
from pathlib import Path

import pytest
from jsonschema.exceptions import SchemaError, ValidationError
from jsonschema.validators import validate
from ruamel.yaml import YAML

from cwl_utils.inputs_schema_gen import cwl_to_jsonschema
from cwl_utils.loghandler import _logger as _cwlutilslogger
from cwl_utils.parser import load_document_by_uri

from .util import get_path

TEST_PARAMS = [
    # Packed Case
    (
        get_path("testdata/revsort-packed.cwl"),
        get_path("testdata/revsort-job.json"),
    ),
    # The number of parameters is a little large, and the definition itself is a straightforward case.
    (
        get_path("testdata/bwa-mem-tool.cwl"),
        get_path("testdata/bwa-mem-job.json"),
    ),
    # The case where CommandInputParameter is shortened (e.g., param: string)
    (
        get_path("testdata/env-tool1.cwl"),
        get_path("testdata/env-job.json"),
    ),
    # Dir
    (
        get_path("testdata/dir.cwl"),
        get_path("testdata/dir-job.yml"),
    ),
    # SecondaryFiles
    (
        get_path("testdata/rename-inputs.cwl"),
        get_path("testdata/rename-inputs.yml"),
    ),
    # Stage array
    (
        get_path("testdata/stage-array.cwl"),
        get_path("testdata/stage-array-job.json"),
    ),
]


@pytest.mark.parametrize("tool_path,inputs_path", TEST_PARAMS)
def test_cwl_inputs_to_jsonschema(tool_path: Path, inputs_path: Path) -> None:
    cwl_obj = load_document_by_uri(tool_path.as_uri())

    _cwlutilslogger.info(f"Generating schema for {tool_path.name}")
    json_schema = cwl_to_jsonschema(cwl_obj)

    _cwlutilslogger.info(
        f"Testing {inputs_path.name} against schema generated for input {tool_path.name}"
    )

    yaml = YAML()

    input_obj = yaml.load(inputs_path)

    try:
        validate(input_obj, json_schema)
    except (ValidationError, SchemaError) as err:
        _cwlutilslogger.error(
            f"Validation failed for {inputs_path.name} "
            f"against schema generated for input {tool_path.name}"
        )
        raise SchemaError(f"{inputs_path.name} failed schema validation") from err


def test_cwl_inputs_to_jsonschema_fails() -> None:
    """Compare tool schema of param 1 against input schema of param 2."""
    tool_path: Path = TEST_PARAMS[0][0]
    inputs_path: Path = TEST_PARAMS[3][1]

    cwl_obj = load_document_by_uri(tool_path.as_uri())

    _cwlutilslogger.info(f"Generating schema for {tool_path.name}")
    json_schema = cwl_to_jsonschema(cwl_obj)

    _cwlutilslogger.info(
        f"Testing {inputs_path.name} against schema generated for input {tool_path.name}"
    )

    yaml = YAML()

    input_obj = yaml.load(inputs_path)

    # We expect this to fail
    with pytest.raises(ValidationError):
        validate(input_obj, json_schema)
