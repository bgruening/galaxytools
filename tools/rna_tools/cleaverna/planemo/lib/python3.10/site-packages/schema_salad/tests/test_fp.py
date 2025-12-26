import pytest

from schema_salad.avro.schema import Names
from schema_salad.exceptions import ValidationException
from schema_salad.schema import load_and_validate, load_schema

from .util import get_data


def test_fp() -> None:
    path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)
    for t in (
        "foreign/foreign_prop1.cwl",
        "foreign/foreign_prop2.cwl",
        "foreign/foreign_prop3.cwl",
        "foreign/foreign_prop4.cwl",
        "foreign/foreign_prop5.cwl",
        "foreign/foreign_prop6.cwl",
        "foreign/foreign_prop7.cwl",
    ):
        path2 = get_data("tests/" + t)
        load_and_validate(
            document_loader,
            avsc_names,
            path2,
            True,
            strict_foreign_properties=False,
        )

    for t in (
        "foreign/foreign_prop1.cwl",
        "foreign/foreign_prop2.cwl",
        "foreign/foreign_prop4.cwl",
        "foreign/foreign_prop5.cwl",
    ):
        path3 = get_data("tests/" + t)
        with pytest.raises(ValidationException):
            try:
                print(t)
                load_and_validate(
                    document_loader,
                    avsc_names,
                    path3,
                    True,
                    strict_foreign_properties=True,
                )
            except ValidationException as e:
                print("\n", e)
                raise
