import json
from typing import Any

import pytest

import schema_salad.metaschema as cg_metaschema
from schema_salad.exceptions import ValidationException
from schema_salad.utils import yaml_no_ts

from .matcher import JsonDiffMatcher
from .util import get_data_uri, get_path


def test_load() -> None:
    doc = {
        "type": "record",
        "fields": [{"name": "hello", "doc": "Hello test case", "type": "string"}],
    }
    rs = cg_metaschema.RecordSchema.fromDoc(
        doc, "http://example.com/", cg_metaschema.LoadingOptions(no_link_check=True)
    )
    assert "record" == rs.type_
    assert rs.fields and "http://example.com/#hello" == rs.fields[0].name
    assert "Hello test case" == rs.fields[0].doc
    assert "string" == rs.fields[0].type_
    assert {
        "type": "record",
        "fields": [
            {
                "name": "http://example.com/#hello",
                "doc": "Hello test case",
                "type": "string",
            }
        ],
    } == rs.save()


def test_err() -> None:
    doc = {"doc": "Hello test case", "type": "string"}
    with pytest.raises(ValidationException):
        cg_metaschema.RecordField.fromDoc(doc, "", cg_metaschema.LoadingOptions())


def test_include() -> None:
    doc = {"name": "hello", "doc": [{"$include": "hello.txt"}], "type": "documentation"}
    path_uri = get_data_uri("tests") + "/_"
    rf = cg_metaschema.Documentation.fromDoc(
        doc,
        "http://example.com/",
        cg_metaschema.LoadingOptions(fileuri=path_uri, no_link_check=True),
    )
    assert "http://example.com/#hello" == rf.name
    assert ["hello world!\n"] == rf.doc
    assert "documentation" == rf.type_
    assert {
        "name": "http://example.com/#hello",
        "doc": ["hello world!\n"],
        "type": "documentation",
    } == rf.save()


def test_import() -> None:
    doc = {"type": "record", "fields": [{"$import": "hellofield.yml"}]}
    lead = get_data_uri("tests")
    rs = cg_metaschema.RecordSchema.fromDoc(
        doc, "http://example.com/", cg_metaschema.LoadingOptions(fileuri=lead + "/_")
    )
    assert "record" == rs.type_
    assert rs.fields and lead + "/hellofield.yml#hello" == rs.fields[0].name
    assert "hello world!\n" == rs.fields[0].doc
    assert "string" == rs.fields[0].type_
    assert {
        "type": "record",
        "fields": [
            {
                "name": lead + "/hellofield.yml#hello",
                "doc": "hello world!\n",
                "type": "string",
            }
        ],
    } == rs.save()


def test_import2() -> None:
    path_uri = get_data_uri("tests/docimp/d1.yml")
    rs = cg_metaschema.load_document(path_uri, "", cg_metaschema.LoadingOptions())
    path2_uri = get_data_uri("tests/docimp/d1.yml")
    assert [
        {
            "doc": [
                "*Hello*",
                "hello 2",
                "*dee dee dee five*",
                "hello 3",
                "hello 4",
                "*dee dee dee five*",
                "hello 5",
            ],
            "type": "documentation",
            "name": path2_uri + "#Semantic_Annotations_for_Linked_Avro_Data",
        }
    ] == [r.save() for r in rs]


def test_err2() -> None:
    doc = {
        "type": "rucord",
        "fields": [{"name": "hello", "doc": "Hello test case", "type": "string"}],
    }
    with pytest.raises(ValidationException):
        cg_metaschema.RecordSchema.fromDoc(doc, "", cg_metaschema.LoadingOptions())


def test_idmap() -> None:
    doc = {
        "type": "record",
        "fields": {"hello": {"doc": "Hello test case", "type": "string"}},
    }
    rs = cg_metaschema.RecordSchema.fromDoc(
        doc, "http://example.com/", cg_metaschema.LoadingOptions(no_link_check=True)
    )
    assert "record" == rs.type_
    assert rs.fields and "http://example.com/#hello" == rs.fields[0].name
    assert "Hello test case" == rs.fields[0].doc
    assert "string" == rs.fields[0].type_
    assert {
        "type": "record",
        "fields": [
            {
                "name": "http://example.com/#hello",
                "doc": "Hello test case",
                "type": "string",
            }
        ],
    } == rs.save()


def test_idmap2() -> None:
    doc = {"type": "record", "fields": {"hello": "string"}}
    rs = cg_metaschema.RecordSchema.fromDoc(
        doc, "http://example.com/", cg_metaschema.LoadingOptions(no_link_check=True)
    )
    assert "record" == rs.type_
    assert rs.fields and "http://example.com/#hello" == rs.fields[0].name
    assert rs.fields[0].doc is None
    assert "string" == rs.fields[0].type_
    assert {
        "type": "record",
        "fields": [{"name": "http://example.com/#hello", "type": "string"}],
    } == rs.save()


def test_load_pt() -> None:
    path_uri = get_data_uri("tests/pt.yml")
    doc = cg_metaschema.load_document(
        path_uri, "", cg_metaschema.LoadingOptions(no_link_check=True)
    )
    assert [
        "https://w3id.org/cwl/salad#null",
        "http://www.w3.org/2001/XMLSchema#boolean",
        "http://www.w3.org/2001/XMLSchema#int",
        "http://www.w3.org/2001/XMLSchema#long",
        "http://www.w3.org/2001/XMLSchema#float",
        "http://www.w3.org/2001/XMLSchema#double",
        "http://www.w3.org/2001/XMLSchema#string",
    ] == doc.symbols


def test_shortname() -> None:
    """Test shortname() function."""
    assert cg_metaschema.shortname("http://example.com/foo") == "foo"
    assert cg_metaschema.shortname("http://example.com/#bar") == "bar"
    assert cg_metaschema.shortname("http://example.com/foo/bar") == "bar"
    assert cg_metaschema.shortname("http://example.com/foo#bar") == "bar"
    assert cg_metaschema.shortname("http://example.com/#foo/bar") == "bar"
    assert cg_metaschema.shortname("http://example.com/foo#bar/baz") == "baz"


@pytest.fixture
def metaschema_pre() -> Any:
    """Prep-parsed schema for testing."""
    with get_path("tests/metaschema-pre.yml").open() as f:
        pre = json.load(f)
    return pre


def test_load_metaschema(metaschema_pre: Any) -> None:
    path_uri = get_data_uri("metaschema/metaschema.yml")
    doc = cg_metaschema.load_document(
        path_uri,
        "",
        cg_metaschema.LoadingOptions(no_link_check=True),
    )
    saved = [d.save(relative_uris=False) for d in doc]
    # with open(get_data("tests/metaschema-pre.yml"), "w") as fh
    #    json.dump(saved, fh, indent=4)
    assert saved == JsonDiffMatcher(metaschema_pre)


def test_load_by_yaml_metaschema(metaschema_pre: Any) -> None:
    path = get_path("metaschema/metaschema.yml")
    with path.open() as path_handle:
        yaml = yaml_no_ts()
        yaml_doc = yaml.load(path_handle)
    doc = cg_metaschema.load_document_by_yaml(
        yaml_doc,
        path.as_uri(),
        None,
    )
    saved = [d.save(relative_uris=False) for d in doc]
    assert saved == JsonDiffMatcher(metaschema_pre)


def test_load_cwlschema() -> None:
    path_uri = get_data_uri("tests/test_schema/CommonWorkflowLanguage.yml")
    doc = cg_metaschema.load_document(
        path_uri,
        "",
        cg_metaschema.LoadingOptions(no_link_check=True),
    )
    path2 = get_path("tests/cwl-pre.yml")
    saved = [d.save(relative_uris=False) for d in doc]
    # with path2.open("w") as f:
    #     json.dump(saved, f, indent=2)
    with path2.open() as f:
        pre = json.load(f)
    assert saved == JsonDiffMatcher(pre)
