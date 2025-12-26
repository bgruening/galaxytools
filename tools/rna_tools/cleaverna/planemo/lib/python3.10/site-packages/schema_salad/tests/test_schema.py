from pathlib import Path
from typing import Any, cast

from schema_salad import schema

from .util import get_data_uri

cwl_file_uri = get_data_uri("tests/test_schema/CommonWorkflowLanguage.yml")


def test_extend_and_specialize_enums(tmp_path: Path) -> None:
    document_loader, _, _, metaschema_loader = schema.load_schema(cwl_file_uri)
    schema_raw_doc = metaschema_loader.fetch(cwl_file_uri)
    schema_doc, _ = metaschema_loader.resolve_all(schema_raw_doc, cwl_file_uri)

    j = schema.extend_and_specialize(cast(list[dict[str, Any]], schema_doc), document_loader)
    CWLType = next((x for x in j if x["name"] == "https://w3id.org/cwl/cwl#CWLType"), None)
    assert CWLType is not None
    symbols = [
        "https://w3id.org/cwl/salad#null",
        "http://www.w3.org/2001/XMLSchema#boolean",
        "http://www.w3.org/2001/XMLSchema#int",
        "http://www.w3.org/2001/XMLSchema#long",
        "http://www.w3.org/2001/XMLSchema#float",
        "http://www.w3.org/2001/XMLSchema#double",
        "http://www.w3.org/2001/XMLSchema#string",
        "https://w3id.org/cwl/cwl#File",
        "https://w3id.org/cwl/cwl#Directory",
    ]
    for symbol in symbols:
        assert symbol in CWLType["symbols"]
