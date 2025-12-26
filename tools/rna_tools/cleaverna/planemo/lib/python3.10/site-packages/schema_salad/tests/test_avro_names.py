"""Avro related tests."""

from schema_salad.avro.schema import Names
from schema_salad.schema import load_schema

from .util import get_data


def test_avro_loading() -> None:
    """Confirm conversion of SALAD style names to avro."""
    path = get_data("tests/test_schema/avro_naming.yml")
    document_loader, avsc_names, schema_metadata, metaschema_loader = load_schema(path)
    assert isinstance(avsc_names, Names)
    assert avsc_names.get_name("com.example.derived_schema.ExtendedThing", None)
