"""Tests of helpful error messages."""

import re

import pytest

import schema_salad
import schema_salad.main
from schema_salad.avro.schema import Names
from schema_salad.exceptions import ValidationException
from schema_salad.ref_resolver import Loader
from schema_salad.schema import load_and_validate, load_schema
from schema_salad.sourceline import cmap

from .util import get_data, get_data_uri


def test_errors() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    for t in (
        "test_schema/test1.cwl",
        "test_schema/test2.cwl",
        "test_schema/test3.cwl",
        "test_schema/test4.cwl",
        "test_schema/test5.cwl",
        "test_schema/test6.cwl",
        "test_schema/test7.cwl",
        "test_schema/test8.cwl",
        "test_schema/test9.cwl",
        "test_schema/test10.cwl",
        "test_schema/test11.cwl",
        "test_schema/test15.cwl",
    ):
        path2 = get_data("tests/" + t)
        with pytest.raises(ValidationException):
            try:
                load_and_validate(
                    document_loader,
                    avsc_names,
                    path2,
                    True,
                )
            except ValidationException as e:
                print("\n", e)
                raise


def test_error_message1() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test1.cwl"
    match = (
        r"""^.+test1\.cwl:2:1: Object\s+'.+test1\.cwl'\s+is\s+not valid """
        + r"""because\s+tried 'Workflow'\s+but
\s+\* missing\s+required\s+field\s+'inputs'
\s+\* missing\s+required\s+field\s+'outputs'
\s+\* missing\s+required\s+field\s+'steps'$"""
    )
    path2 = get_data("tests/" + t)
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, path2, True)


def test_error_message2() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test2.cwl"
    match = (
        r"""
^.+test2\.cwl:2:1: Field """
        r"""'class'\s+contains\s+undefined\s+reference\s+to\s+'file://.+/schema_salad/tests/test_schema/xWorkflow'$"""[
            1:
        ]
    )
    path2 = get_data("tests/" + t)
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, path2, True)


def test_error_message3() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test3.cwl"
    match = r"""
^.+test3\.cwl:6:1: checking field\s+'outputs'
.+test3\.cwl:7:3:   checking object\s+'.+test3\.cwl#bar'
\s+Field 'type'\s+references\s+unknown\s+identifier\s+'xstring',\s+tried\s+file://.+/tests/test_schema/test3\.cwl#xstring$"""[  # noqa: B950
        1:
    ]
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message4() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test4.cwl"
    match = r"""
^.+test4\.cwl:6:1: checking field\s+'outputs'
.+test4\.cwl:7:3:   checking object\s+'.+test4\.cwl#bar'
\s+'type'\s+field\s+is\s+int,\s+expected\s+string,\s+list,\s+or\s+a\s+dict.$"""[
        1:
    ]
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message5() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test5.cwl"
    match = r"""
^.+test5\.cwl:2:1: Object\s+'.+test5\.cwl'\s+is\s+not valid because
\s+tried 'Workflow'\s+but
.+test5\.cwl:8:1:\s+the 'steps'\s+field\s+is\s+not\s+valid\s+because
\s+tried array\s+of\s+<WorkflowStep>\s+but
.+test5\.cwl:8:9:\s+item is\s+invalid\s+because
\s+is not a\s+dict.\s+Expected\s+a\s+WorkflowStep\s+object.$"""[
        1:
    ]
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message7() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test7.cwl"
    match = (
        r"""^.+test7\.cwl:2:1:\s+Object\s+'.+test7\.cwl'\s+is\s+not\s+valid\s+because
\s+tried\s+'Workflow'\s+but
.+test7\.cwl:8:1:\s+the 'steps'\s+field\s+is\s+not\s+valid\s+because
\s+tried\s+array\s+of\s+<WorkflowStep>\s+but
.+test7\.cwl:9:3:\s+item is\s+invalid\s+because
\s+\* missing\s+required\s+field\s+'run'
.+test7\.cwl:10:5:\s+\* invalid\s+field\s+'scatter_method',\s+expected\s+one\s+"""
        + r"of:\s+'id',\s+'in',\s+'out',\s+'requirements',\s+'hints',\s+"
        + r"'label',\s+'doc',\s+'run',\s+'scatter',\s+'scatterMethod'$"
    )
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message8() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test8.cwl"
    match = (
        r"""
^.+test8\.cwl:8:1:\s+checking field\s+'steps'
.+test8\.cwl:9:3:\s+checking object\s+'.+test8\.cwl#step1'
.+test8\.cwl:10:5:\s+"""
        r"""Field\s+'scatterMethod'\s+contains\s+undefined\s+reference\s+to\s+'file:///.+/tests/test_schema/abc'$"""[
            1:
        ]
    )
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message9() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test9.cwl"
    match = (
        r"""^.+test9\.cwl:8:1:\s+checking field\s+'steps'
.+test9\.cwl:9:3:\s+checking object\s+'.+test9\.cwl#step1'
.+test9\.cwl:10:5:\s+'scatterMethod'\s+field\s+is\s+"""
        + r"""int,\s+expected\s+string,\s+list,\s+or a\s+dict.$"""
    )
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message10() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test10.cwl"
    match = r"""
^.+test10\.cwl:2:1:\s+Object\s+'.+test10\.cwl'\s+is\s+not\s+valid\s+because
\s+tried 'Workflow'\s+but
.+test10\.cwl:8:1:\s+the 'steps'\s+field\s+is\s+not\s+valid\s+because
\s+tried array\s+of\s+<WorkflowStep>\s+but
.+test10\.cwl:9:3:\s+item is\s+invalid\s+because
\s+\* missing\s+required\s+field\s+'run'
.+test10\.cwl:10:5:\s+\* the\s+'scatterMethod'\s+field\s+is\s+not\s+valid\s+because
\s+value\s+is\s+a\s+CommentedSeq,\s+expected\s+null\s+or\s+ScatterMethod$"""[
        1:
    ]
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message11() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test11.cwl"
    match = (
        r"""
^.+test11\.cwl:8:1:\s+checking field\s+'steps'
.+test11\.cwl:9:3:\s+checking object\s+'.+test11\.cwl#step1'
.+test11\.cwl:10:5:\s+"""
        r"""Field 'run'\s+contains\s+undefined\s+reference\s+to\s+'file://.+/tests/test_schema/blub\.cwl'$"""[
            1:
        ]
    )
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


def test_error_message15() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    t = "test_schema/test15.cwl"
    match = (
        r"""^.+test15\.cwl:3:1:\s+Object\s+'.+test15\.cwl'\s+is\s+not\s+valid\s+because
\s+tried\s+'CommandLineTool'\s+but
.+test15\.cwl:6:1:\s+the\s+'inputs'\s+field\s+is\s+not\s+valid\s+because
.+test15\.cwl:7:3:\s+item\s+is\s+invalid\s+because
.+test15\.cwl:9:5:\s+the\s+'inputBinding'\s+field\s+is\s+not\s+valid\s+because
.+tried\s+CommandLineBinding\s+but
.+test15\.cwl:11:7:             \*\s+invalid\s+field\s+'invalid_field',\s+expected\s+"""
        + r"""one\s+of:\s+'loadContents',\s+'position',\s+'prefix',\s+'separate',"""
        + r"""\s+'itemSeparator',\s+'valueFrom',\s+'shellQuote'
.+test15\.cwl:12:7:             \*\s+invalid\s+field\s+'another_invalid_field',"""
        + r"""\s+expected\s+one\s+of:\s+'loadContents',\s+'position',\s+'prefix',"""
        + r"""\s+'separate',\s+'itemSeparator',\s+'valueFrom',\s+'shellQuote'$"""
    )
    with pytest.raises(ValidationException, match=match):
        load_and_validate(document_loader, avsc_names, get_data("tests/" + t), True)


@pytest.mark.skip(
    "See https://github.com/common-workflow-language/common-workflow-language/issues/734"  # noqa: B950
)
def test_errors_previously_defined_dict_key() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)

    for t in (
        "test_schema/test12.cwl",
        "test_schema/test13.cwl",
        "test_schema/test14.cwl",
    ):
        with pytest.raises(ValidationException):
            try:
                load_and_validate(
                    document_loader,
                    avsc_names,
                    get_data("tests/" + t),
                    True,
                )
            except ValidationException as e:
                print("\n", e)
                raise


def test_bad_schema() -> None:
    path = get_data("tests/bad_schema.yml")
    assert 1 == schema_salad.main.main(argsl=[path])
    path2 = get_data("tests/bad_schema.yml")
    assert 1 == schema_salad.main.main(argsl=["--print-avro", path2])


def test_bad_schema2() -> None:
    path = get_data("tests/bad_schema2.yml")
    assert 1 == schema_salad.main.main(argsl=[path])


def test_namespaces_type() -> None:
    """Confirm helpful error message when $namespaces is the wrong type."""
    with pytest.raises(
        ValidationException,
        match=re.escape("test:1:1: $namespaces must be a dictionary"),
    ):
        ldr, _, _, _ = schema_salad.schema.load_schema(
            cmap(
                {
                    "$base": "Y",
                    "name": "X",
                    "$namespaces": ["http://example.com/foo#"],
                    "$graph": [
                        {
                            "name": "ExampleType",
                            "type": "enum",
                            "symbols": ["asym", "bsym"],
                        }
                    ],
                }
            )
        )


def test_namespaces_undeclared(caplog: pytest.LogCaptureFixture) -> None:
    """Confirm warning message a namespace is used but not declared."""

    ldr, _, _, _ = schema_salad.schema.load_schema(
        cmap(
            {
                "$base": "Y",
                "name": "X",
                "$graph": [
                    {
                        "name": "namesp:ExampleType",
                        "type": "enum",
                        "symbols": ["asym", "bsym"],
                    }
                ],
            }
        )
    )

    match = (
        r""".*URI prefix 'namesp' of 'namesp:ExampleType' not recognized, """
        r"""are you missing a \$namespaces section?.*"""
    )

    assert re.match(match, caplog.text)


def test_not_a_namespace1(caplog: pytest.LogCaptureFixture) -> None:
    """Confirm no warning when relative id contains a colon but prefix doesn't look like a namespace."""

    ldr, _, _, _ = schema_salad.schema.load_schema(
        cmap(
            {
                "$base": "Y",
                "name": "X",
                "$graph": [
                    {
                        "name": "foo/colon:ExampleType",
                        "type": "enum",
                        "symbols": ["asym", "bsym"],
                    }
                ],
            }
        )
    )

    assert caplog.text == ""


def test_not_a_namespace2(caplog: pytest.LogCaptureFixture) -> None:
    """Confirm no warning when relative id contains a colon but prefix doesn't look like a namespace."""

    ldr, _, _, _ = schema_salad.schema.load_schema(
        cmap(
            {
                "$base": "Y",
                "name": "X",
                "$graph": [
                    {
                        "name": "foo#colon:ExampleType",
                        "type": "enum",
                        "symbols": ["asym", "bsym"],
                    }
                ],
            }
        )
    )

    assert caplog.text == ""


def test_not_a_namespace3(caplog: pytest.LogCaptureFixture) -> None:
    """Confirm no warning when relative id starts with a colon."""

    ldr, _, _, _ = schema_salad.schema.load_schema(
        cmap(
            {
                "$base": "Y",
                "name": "X",
                "$graph": [
                    {
                        "name": ":ExampleType",
                        "type": "enum",
                        "symbols": ["asym", "bsym"],
                    }
                ],
            }
        )
    )

    assert caplog.text == ""


def test_schemas_type() -> None:
    """Confirm helpful error message when $schemas is the wrong type."""
    path_uri = get_data_uri("tests/EDAM.owl")
    with pytest.raises(
        ValidationException,
        match=re.escape("test:1:1: $schemas must be a string or a list of string"),
    ):
        ra, _ = Loader({}).resolve_all(
            cmap(
                {
                    "$schemas": {"wrong": path_uri},
                    "edam:has_format": "edam:format_1915",
                }
            ),
            "",
        )
