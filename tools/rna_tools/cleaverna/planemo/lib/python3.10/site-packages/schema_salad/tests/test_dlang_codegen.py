"""Test D code generation."""

import os
from pathlib import Path
from typing import Any, cast

from schema_salad import codegen
from schema_salad.avro.schema import Names
from schema_salad.schema import load_schema

from .util import cwl_file_uri


def test_cwl_dlang_gen(tmp_path: Path) -> None:
    """End to end test of D generator using the CWL v1.0 schema."""
    src_target = tmp_path / "cwl_v1_0.d"
    dlang_codegen(cwl_file_uri, src_target)
    assert os.path.exists(src_target)


def dlang_codegen(
    file_uri: str,
    target: Path,
) -> None:
    """Help using the D code generation function."""
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(file_uri)
    assert isinstance(avsc_names, Names)
    schema_raw_doc = metaschema_loader.fetch(file_uri)
    schema_doc, schema_metadata = metaschema_loader.resolve_all(schema_raw_doc, file_uri)
    codegen.codegen(
        "dlang",
        cast(list[dict[str, Any]], schema_doc),
        schema_metadata,
        document_loader,
        target=str(target),
    )
