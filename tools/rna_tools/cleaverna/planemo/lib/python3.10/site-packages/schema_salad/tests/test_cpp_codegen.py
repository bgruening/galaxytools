"""Test C++ code generation."""

import filecmp
import os
from pathlib import Path
from typing import Any, Optional, cast

import pytest

from schema_salad import codegen
from schema_salad.avro.schema import Names
from schema_salad.schema import load_schema

from .util import cwl_file_uri, get_data_uri, get_path


def test_cwl_cpp_gen(tmp_path: Path) -> None:
    """End to end test of C++ generator using the CWL v1.0 schema."""
    src_target = tmp_path / "cwl_v1_0.h"
    cpp_codegen(cwl_file_uri, src_target)
    assert os.path.exists(src_target)


@pytest.mark.parametrize(
    "filename",
    [
        "01_single_record.yml",
        "02_two_records.yml",
        "03_simple_inheritance.yml",
        "04_abstract_inheritance.yml",
        "05_specialization.yml",
    ],
)
def test_cwl_cpp_generations(tmp_path: Path, filename: str) -> None:
    """End to end test of C++ generator using small scenarios."""

    file = get_path(f"cpp_tests/{filename}")

    # file with generated cpp output
    src_target = tmp_path / "test.h"
    # file with expected cpp output
    expected = file.with_suffix(".h")

    cpp_codegen("file://" + os.fspath(file), src_target)

    assert os.path.isfile(expected)
    assert os.path.isfile(src_target)
    assert filecmp.cmp(expected, src_target, shallow=False)


def test_cwl_cpp_generations_with_spdx(tmp_path: Path) -> None:
    """End to end test of C++ generator checking for SPDX headers"""

    src_target = tmp_path / "test.h"

    input_file_uri = get_data_uri("cpp_tests/01_single_record.yml")

    """Generating different combinations of license headers"""
    """Generate License Identifier"""
    cpp_codegen(input_file_uri, src_target, spdx_license_identifier="Apache-2.0")
    lines = open(src_target).readlines()[0:2]
    assert lines[0] == "// SPDX-License-Identifier: Apache-2.0\n"
    assert lines[1] == "#pragma once\n"

    """Generate single CopyrightText"""
    cpp_codegen(
        input_file_uri,
        src_target,
        spdx_copyright_text=["Copyright 2016 Some People <email@example.com>"],
    )
    lines = open(src_target).readlines()[0:2]
    assert lines[0] == "// SPDX-FileCopyrightText: Copyright 2016 Some People <email@example.com>\n"
    assert lines[1] == "#pragma once\n"

    """Generate two CopyrightText entries"""
    cpp_codegen(
        input_file_uri,
        src_target,
        spdx_copyright_text=[
            "Copyright 2016 Person A <person_a@example.com>",
            "Copyright 2017 Person B <person_b@example.com>",
        ],
    )
    lines = open(src_target).readlines()[0:3]
    assert lines[0] == "// SPDX-FileCopyrightText: Copyright 2016 Person A <person_a@example.com>\n"
    assert lines[1] == "// SPDX-FileCopyrightText: Copyright 2017 Person B <person_b@example.com>\n"
    assert lines[2] == "#pragma once\n"

    """Generate CopyrightText and License Identifier"""
    cpp_codegen(
        input_file_uri,
        src_target,
        spdx_license_identifier="Apache-2.0",
        spdx_copyright_text=[
            "Copyright 2016 Person A <person_a@example.com>",
            "Copyright 2017 Person B <person_b@example.com>",
        ],
    )
    lines = open(src_target).readlines()[0:4]

    assert lines[0] == "// SPDX-FileCopyrightText: Copyright 2016 Person A <person_a@example.com>\n"
    assert lines[1] == "// SPDX-FileCopyrightText: Copyright 2017 Person B <person_b@example.com>\n"
    assert lines[2] == "// SPDX-License-Identifier: Apache-2.0\n"
    assert lines[3] == "#pragma once\n"


def cpp_codegen(
    file_uri: str,
    target: Path,
    spdx_copyright_text: Optional[list[str]] = None,
    spdx_license_identifier: Optional[str] = None,
) -> None:
    """Help using the C++ code generation function."""
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(file_uri)
    assert isinstance(avsc_names, Names)
    schema_raw_doc = metaschema_loader.fetch(file_uri)
    schema_doc, schema_metadata = metaschema_loader.resolve_all(schema_raw_doc, file_uri)
    codegen.codegen(
        "cpp",
        cast(list[dict[str, Any]], schema_doc),
        schema_metadata,
        document_loader,
        target=str(target),
        spdx_copyright_text=spdx_copyright_text,
        spdx_license_identifier=spdx_license_identifier,
    )
