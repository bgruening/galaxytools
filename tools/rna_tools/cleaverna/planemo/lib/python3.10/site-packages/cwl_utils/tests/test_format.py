# SPDX-License-Identifier: Apache-2.0
"""Tests of cwl_utils.file_formats."""

import xml.sax
from pathlib import Path
from typing import Optional

import requests
from pytest import raises
from rdflib import Graph
from rdflib.compare import to_isomorphic
from rdflib.plugins.parsers.notation3 import BadSyntax
from schema_salad.exceptions import ValidationException
from schema_salad.fetcher import DefaultFetcher

from cwl_utils.file_formats import check_format
from cwl_utils.parser import load_document_by_uri
from cwl_utils.types import CWLObjectType

from .util import get_path


def _create_file(format_: Optional[str] = None) -> CWLObjectType:
    obj: CWLObjectType = {
        "class": "File",
        "basename": "example.txt",
        "size": 23,
        "contents": "hoopla",
        "nameroot": "example",
        "nameext": "txt",
    }
    if format_:
        obj["format"] = format_
    return obj


def _load_format(fetchpath: Path) -> Graph:
    fetcher = DefaultFetcher({}, requests.Session())
    fetchurl = fetchpath.as_uri()
    content = fetcher.fetch_text(fetchurl)
    graph = Graph()
    for fmt in ["xml", "turtle", "rdfa"]:
        try:
            graph.parse(data=content, format=fmt, publicID=str(fetchurl))
            break
        except (xml.sax.SAXParseException, TypeError, BadSyntax):
            pass
    return graph


EDAM = _load_format(get_path("testdata/EDAM_subset.owl"))
GX = _load_format(get_path("testdata/gx_edam.ttl"))


def test_check_format() -> None:
    """Exact format equivalence test, with ontology."""
    check_format(
        actual_file=_create_file(format_="http://edamontology.org/format_2330"),
        input_formats="http://edamontology.org/format_2330",
        ontology=EDAM,
    )


def test_check_format_subformat() -> None:
    """Test of check_format with a subformat."""
    check_format(
        actual_file=_create_file(format_="http://edamontology.org/format_1929"),
        input_formats="http://edamontology.org/format_2330",
        ontology=EDAM,
    )


def test_check_format_equiv() -> None:
    """Test of check_format with an equivalent format."""
    check_format(
        actual_file=_create_file(format_="http://edamontology.org/format_1929"),
        input_formats="http://galaxyproject.org/formats/fasta",
        ontology=EDAM + GX,
    )


def test_check_format_equiv2() -> None:
    """Test of check_format with an equivalent format, in the reverse."""
    check_format(
        actual_file=_create_file(format_="http://galaxyproject.org/formats/fasta"),
        input_formats="http://edamontology.org/format_1929",
        ontology=EDAM + GX,
    )


def test_check_format_wrong_format() -> None:
    """Test of check_format with a non-match format with an ontology."""
    with raises(ValidationException, match=r"File has an incompatible format: .*"):
        check_format(
            actual_file=_create_file(format_="http://edamontology.org/format_1929"),
            input_formats="http://edamontology.org/format_2334",
            ontology=EDAM,
        )


def test_check_format_wrong_format_no_ontology() -> None:
    """Test of check_format with a non-match format."""
    with raises(ValidationException, match=r"File has an incompatible format: .*"):
        check_format(
            actual_file=_create_file(format_="http://edamontology.org/format_1929"),
            input_formats="http://edamontology.org/format_2334",
            ontology=None,
        )


def test_check_format_no_format() -> None:
    """Confirm that a missing format produces the expected exception."""
    with raises(ValidationException, match=r"File has no 'format' defined: .*"):
        check_format(
            actual_file=_create_file(),
            input_formats="http://edamontology.org/format_2330",
            ontology=EDAM,
        )


def test_check_format_missing_file() -> None:
    """Confirm that a missing file produces no error."""
    check_format(
        actual_file=[{}],
        input_formats="http://edamontology.org/format_2330",
        ontology=EDAM,
    )


def test_check_format_no_ontology() -> None:
    """Confirm that precisely matching formats without an ontology still match."""
    check_format(
        actual_file=_create_file(format_="http://edamontology.org/format_2330"),
        input_formats="http://edamontology.org/format_2330",
        ontology=Graph(),
    )


def test_loading_options_graph_property_v1_0() -> None:
    """Test that RDFLib Graph representations of $schema properties are correctly loaded, CWL v1.0."""
    uri = get_path("testdata/formattest2_v1_0.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert to_isomorphic(cwl_obj.loadingOptions.graph) == to_isomorphic(EDAM)


def test_loading_options_graph_property_v1_1() -> None:
    """Test that RDFLib Graph representations of $schema properties are correctly loaded, CWL v1.1."""
    uri = get_path("testdata/formattest2_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert to_isomorphic(cwl_obj.loadingOptions.graph) == to_isomorphic(EDAM)


def test_loading_options_graph_property_v1_2() -> None:
    """Test that RDFLib Graph representations of $schema properties are correctly loaded, CWL v1.2."""
    uri = get_path("testdata/formattest2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert to_isomorphic(cwl_obj.loadingOptions.graph) == to_isomorphic(EDAM)


def test_loading_options_missing_graph_v1_0() -> None:
    """Affirm that v1.0 documents without $schema still produce an empty graph property."""
    uri = get_path("testdata/workflow_input_format_expr.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert to_isomorphic(cwl_obj.loadingOptions.graph) == to_isomorphic(Graph())


def test_loading_options_missing_graph_v1_1() -> None:
    """Affirm that v1.1 documents without $schema still produce an empty graph property."""
    uri = get_path("testdata/workflow_input_format_expr_v1_1.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert to_isomorphic(cwl_obj.loadingOptions.graph) == to_isomorphic(Graph())


def test_loading_options_missing_graph_v1_2() -> None:
    """Affirm that v1.2 documents without $schema still produce an empty graph property."""
    uri = get_path("testdata/workflow_input_format_expr_v1_2.cwl").as_uri()
    cwl_obj = load_document_by_uri(uri)
    assert to_isomorphic(cwl_obj.loadingOptions.graph) == to_isomorphic(Graph())


def test_loading_format_without_schema_v1_0() -> None:
    """Test that format fields without accompanying schemas are tolerated, CWL v1.0."""
    uri = get_path("testdata/revsort-run-1-packed.cwl").as_uri()
    load_document_by_uri(uri)
