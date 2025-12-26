import inspect
import os
from pathlib import Path
from typing import Any, Optional, cast

from rdflib import Graph
from rdflib.compare import to_isomorphic
from requests import Session

import schema_salad.metaschema as cg_metaschema
from schema_salad import codegen
from schema_salad.avro.schema import Names
from schema_salad.fetcher import DefaultFetcher
from schema_salad.python_codegen import PythonCodeGen
from schema_salad.python_codegen_support import LoadingOptions
from schema_salad.schema import load_schema

from .util import basket_file_uri, cwl_file_uri, get_data, get_path, metaschema_file_uri


def test_safe_identifiers() -> None:
    """Affirm correct construction of identifiers safe for Python."""
    assert PythonCodeGen.safe_name("5.8s_pattern") == "_5_8s_pattern"


def test_cwl_gen(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(cwl_file_uri, src_target)
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert "class Workflow(Process)" in f.read()


def test_meta_schema_gen(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(metaschema_file_uri, src_target)
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert "class RecordSchema(Saveable):" in f.read()


def test_meta_schema_gen_up_to_date(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(metaschema_file_uri, src_target)
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert f.read() == inspect.getsource(cg_metaschema)


def test_meta_schema_gen_no_base(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(basket_file_uri, src_target)
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert "class Basket" in f.read()


def python_codegen(
    file_uri: str,
    target: Path,
    parser_info: Optional[str] = None,
    package: Optional[str] = None,
) -> None:
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(file_uri)
    assert isinstance(avsc_names, Names)
    schema_raw_doc = metaschema_loader.fetch(file_uri)
    schema_doc, schema_metadata = metaschema_loader.resolve_all(schema_raw_doc, file_uri)
    codegen.codegen(
        "python",
        cast(list[dict[str, Any]], schema_doc),
        schema_metadata,
        document_loader,
        target=str(target),
        parser_info=parser_info,
        package=package,
    )


def test_default_parser_info(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(metaschema_file_uri, src_target)
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert 'def parser_info() -> str:\n    return "org.w3id.cwl.salad"' in f.read()


def test_parser_info(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(metaschema_file_uri, src_target, parser_info="cwl")
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert 'def parser_info() -> str:\n    return "cwl"' in f.read()


def test_use_of_package_for_parser_info(tmp_path: Path) -> None:
    src_target = tmp_path / "src.py"
    python_codegen(metaschema_file_uri, src_target, package="cwl")
    assert os.path.exists(src_target)
    with open(src_target) as f:
        assert 'def parser_info() -> str:\n    return "cwl"' in f.read()


def test_graph_property() -> None:
    """Test the RDFLib Graph representation of the `$schemas` directive."""
    schema = get_path("tests/EDAM.owl")
    fetcher = DefaultFetcher({}, Session())
    fetchurl = schema.as_uri()
    content = fetcher.fetch_text(fetchurl)
    graph = Graph()
    graph.parse(data=content, format="xml", publicID=fetchurl)
    loading_options = LoadingOptions(schemas=[str(schema)])
    assert to_isomorphic(graph) == to_isomorphic(loading_options.graph)


def test_graph_property_cache() -> None:
    """Test that LoadingOptions properly cache the `$schemas` RDFLib Graph representations."""
    schema = get_data("tests/EDAM.owl")
    loading_options = LoadingOptions(schemas=[schema])
    graph1 = loading_options.graph
    graph2 = loading_options.graph
    assert graph1 == graph2


def test_graph_property_empty_schema() -> None:
    """Test that an empty RDFLib Graph is returned when not `$schemas` directive is present."""
    loading_options = LoadingOptions()
    assert to_isomorphic(loading_options.graph) == to_isomorphic(Graph())
