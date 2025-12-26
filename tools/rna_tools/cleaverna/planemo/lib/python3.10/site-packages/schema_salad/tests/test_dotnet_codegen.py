import shutil
from pathlib import Path
from typing import Any, Optional, cast

from schema_salad import codegen
from schema_salad.schema import load_schema

from .util import cwl_file_uri, get_data_uri, get_path, metaschema_file_uri


def test_cwl_gen(tmp_path: Path) -> None:
    topmed_example_path = get_path("tests/test_real_cwl/topmed/topmed_variant_calling_pipeline.cwl")
    target_dir = tmp_path / "target"
    examples_dir = tmp_path / "examples"

    target_dir.mkdir()
    examples_dir.mkdir()
    shutil.copyfile(topmed_example_path, examples_dir / "valid_topmed.cwl")

    dotnet_codegen(cwl_file_uri, target_dir, examples=examples_dir)
    solution_path = target_dir / "Solution.sln"
    assert solution_path.exists()
    tests_dir = target_dir / "Test" / "src"
    assert tests_dir.exists()
    with open(tests_dir / "ExampleTests.cs") as f:
        assert "topmed" in f.read()


def test_meta_schema_gen(tmp_path: Path) -> None:
    target_dir = tmp_path / "target"
    target_dir.mkdir()
    dotnet_codegen(metaschema_file_uri, target_dir)
    solution_path = target_dir / "Solution.sln"
    assert solution_path.exists()
    src_dir = target_dir / "DotnetTest" / "src"
    assert src_dir.exists()
    record_schema_dir = src_dir / "RecordSchema.cs"
    assert record_schema_dir.exists()
    with open(record_schema_dir) as f:
        assert "public class RecordSchema : IRecordSchema, " "ISaveable\n{\n" in f.read()


def test_class_field(tmp_path: Path) -> None:
    schema_path = get_data_uri("tests/class_field_test.yml")
    assert schema_path
    target_dir = tmp_path / "target"

    target_dir.mkdir()
    dotnet_codegen(schema_path, target_dir)

    solution_path = target_dir / "Solution.sln"
    assert solution_path.exists()

    tests_dir = target_dir / "Test"
    assert tests_dir.exists()

    with open(target_dir / "DotnetTest" / "src" / "ClassFieldString.cs") as f:
        assert (
            'public ClassFieldString(string class_ = "ClassFieldString",'
            + " LoadingOptions? loadingOptions = null,"
            + " Dictionary<object, object>? extensionFields = null)\n    {\n"
            in f.read()
        )
    with open(target_dir / "DotnetTest" / "src" / "ClassFieldEnum.cs") as f:
        assert (
            "public ClassFieldEnum(ClassFieldEnum_class? class_ = null,"
            + " LoadingOptions? loadingOptions = null,"
            + " Dictionary<object, object>? extensionFields = null)\n    {\n"
            in f.read()
        )


def dotnet_codegen(file_uri: str, target: Path, examples: Optional[Path] = None) -> None:
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(file_uri)
    schema_raw_doc = metaschema_loader.fetch(file_uri)
    schema_doc, schema_metadata = metaschema_loader.resolve_all(schema_raw_doc, file_uri)
    codegen.codegen(
        "dotnet",
        cast(list[dict[str, Any]], schema_doc),
        schema_metadata,
        document_loader,
        package="DotnetTest",
        target=str(target),
        examples=str(examples) if examples else None,
    )
