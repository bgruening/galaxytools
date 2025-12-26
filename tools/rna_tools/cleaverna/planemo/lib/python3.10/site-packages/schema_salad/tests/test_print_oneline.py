import re
from os.path import normpath

import pytest

from schema_salad.avro.schema import Names
from schema_salad.exceptions import ValidationException, to_one_line_messages
from schema_salad.schema import load_and_validate, load_schema

from .util import get_data


def test_print_oneline() -> None:
    # Issue #135
    path1 = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path1)
    assert isinstance(avsc_names, Names)

    src = "test15.cwl"
    path2 = get_data("tests/test_schema/" + src)
    with pytest.raises(ValidationException):
        try:
            load_and_validate(
                document_loader,
                avsc_names,
                path2,
                True,
            )
        except ValidationException as e:
            msgs = to_one_line_messages(e).splitlines()
            assert len(msgs) == 2
            assert msgs[0].endswith(
                src + ":11:7: invalid field 'invalid_field', expected one of: "
                "'loadContents', 'position', 'prefix', 'separate', "
                "'itemSeparator', 'valueFrom', 'shellQuote'"
            )
            assert msgs[1].endswith(
                src + ":12:7: invalid field 'another_invalid_field', expected one of: "
                "'loadContents', 'position', 'prefix', 'separate', 'itemSeparator', "
                "'valueFrom', 'shellQuote'"
            )
            print("\n", e)
            raise


def test_print_oneline_for_invalid_yaml() -> None:
    # Issue #137
    path1 = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path1)
    assert isinstance(avsc_names, Names)

    src = "test16.cwl"
    path2 = get_data("tests/test_schema/" + src)
    with pytest.raises(ValidationException):
        try:
            load_and_validate(
                document_loader,
                avsc_names,
                path2,
                True,
            )
        except ValidationException as e:
            msg = to_one_line_messages(e)
            print("\n", e)
            assert msg.endswith(src + ":11:1: could not find expected ':'")
            raise


def test_print_oneline_for_errors_in_the_same_line() -> None:
    # Issue #136
    path1 = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path1)
    assert isinstance(avsc_names, Names)

    src = "test17.cwl"
    path2 = get_data("tests/test_schema/" + src)
    with pytest.raises(ValidationException):
        try:
            load_and_validate(
                document_loader,
                avsc_names,
                path2,
                True,
            )
        except ValidationException as e:
            msgs = to_one_line_messages(e).splitlines()
            assert len(msgs) == 2, msgs
            assert msgs[0].endswith(src + ":14:5: missing required field 'id'")
            assert msgs[1].endswith(
                src + ":14:5: invalid field 'aa', expected one of: 'label', "
                "'secondaryFiles', 'streamable', 'doc', 'id', "
                "'outputBinding', 'format', 'type'"
            )
            print("\n", e)
            raise


def test_print_oneline_for_errors_in_resolve_ref() -> None:
    # Issue #141
    path1 = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path1)
    assert isinstance(avsc_names, Names)

    src = "test18.cwl"
    path2 = get_data("tests/test_schema/" + src)
    fullpath = normpath(path2)
    with pytest.raises(ValidationException):
        try:
            load_and_validate(document_loader, avsc_names, str(fullpath), True)
        except ValidationException as e:
            msg = to_one_line_messages(e)
            # convert Windows path to Posix path
            if "\\" in fullpath:
                fullpath = "/" + fullpath.lower().replace("\\", "/")
            print("\n", e)
            assert msg.endswith(
                src + ":14:5: Field 'type' references unknown identifier "
                "'Filea', tried file://{}#Filea".format(fullpath)
            )
            raise


def test_for_invalid_yaml1() -> None:
    # Issue 143
    path1 = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path1)
    assert isinstance(avsc_names, Names)

    src = "test16.cwl"
    path2 = get_data("tests/test_schema/" + src)
    with pytest.raises(ValidationException):
        try:
            load_and_validate(
                document_loader,
                avsc_names,
                path2,
                True,
            )
        except ValidationException as e:
            msg = str(e)
            print("\n", e)
            assert re.search(src + r":10:7: while\s+scanning\s+a\s+simple\s+key", msg, re.M)
            assert re.search(src + r":11:1:   could\s+not\s+find\s+expected ':'$", msg, re.M)
            raise


def test_for_invalid_yaml2() -> None:
    # Issue 143
    path1 = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path1)
    assert isinstance(avsc_names, Names)

    src = "test19.cwl"
    path2 = get_data("tests/test_schema/" + src)
    with pytest.raises(ValidationException):
        try:
            load_and_validate(
                document_loader,
                avsc_names,
                path2,
                True,
            )
        except ValidationException as e:
            msg = str(e)
            assert (
                msg.endswith("expected <block end>, but found ':'")
                or msg.endswith("expected <block end>, but found u':'")
                or re.search(r"mapping\s+with\s+implicit\s+null\s+key$", msg, re.M)
            )
            print("\n", e)
            raise
