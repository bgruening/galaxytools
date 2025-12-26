"""
Tests to ensure that mypyc compiled classes are still pickleable.

See https://mypyc.readthedocs.io/en/latest/differences_from_python.html#pickling-and-copying-objects
"""

import pickle
from pathlib import Path

from schema_salad import ref_resolver, schema
from schema_salad.avro.schema import Names, RecordSchema

from .util import get_data_uri


def test_recordschema_pickle() -> None:
    """Targeted test of pickling a RecordSchema."""
    s = RecordSchema("one", None, [], Names())
    print(s)
    d = pickle.dumps(s)
    print(pickle.loads(d))


def test_loader_pickle() -> None:
    """Pickle a Loader."""
    loader = ref_resolver.Loader({})
    print(loader)
    d = pickle.dumps(loader)
    print(pickle.loads(d))


def test_extend_and_specialize_enums(tmp_path: Path) -> None:
    cwl_file_uri = get_data_uri("tests/test_schema/CommonWorkflowLanguage.yml")
    _, avsc_names, _, _ = schema.load_schema(cwl_file_uri)
    print(avsc_names)
    print(pickle.loads(pickle.dumps(avsc_names)))
