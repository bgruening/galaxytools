import shutil
from pathlib import Path
from typing import Any, Optional, cast

from schema_salad import codegen
from schema_salad.schema import load_schema

from .util import cwl_file_uri, get_path, metaschema_file_uri


def test_cwl_gen(tmp_path: Path) -> None:
    topmed_example_path = get_path("tests/test_real_cwl/topmed/topmed_variant_calling_pipeline.cwl")
    target_dir = tmp_path / "target"
    examples_dir = tmp_path / "examples"

    target_dir.mkdir()
    examples_dir.mkdir()
    shutil.copyfile(topmed_example_path, examples_dir / "valid_topmed.cwl")

    java_codegen(cwl_file_uri, target_dir, examples=examples_dir)
    pom_xml_path = target_dir / "pom.xml"
    assert pom_xml_path.exists()
    tests_dir = target_dir / "src" / "test" / "java" / "org" / "w3id" / "cwl" / "cwl" / "utils"
    assert tests_dir.exists()
    with open(tests_dir / "ExamplesTest.java") as f:
        assert "topmed" in f.read()


def test_meta_schema_gen(tmp_path: Path) -> None:
    target_dir = tmp_path / "target"
    target_dir.mkdir()
    java_codegen(metaschema_file_uri, target_dir)
    pom_xml_path = target_dir / "pom.xml"
    assert pom_xml_path.exists()
    src_dir = target_dir / "src" / "main" / "java" / "org" / "w3id" / "cwl" / "salad"
    assert src_dir.exists()


def java_codegen(file_uri: str, target: Path, examples: Optional[Path] = None) -> None:
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(file_uri)
    schema_raw_doc = metaschema_loader.fetch(file_uri)
    schema_doc, schema_metadata = metaschema_loader.resolve_all(schema_raw_doc, file_uri)
    codegen.codegen(
        "java",
        cast(list[dict[str, Any]], schema_doc),
        schema_metadata,
        document_loader,
        target=str(target),
        examples=str(examples) if examples else None,
    )
