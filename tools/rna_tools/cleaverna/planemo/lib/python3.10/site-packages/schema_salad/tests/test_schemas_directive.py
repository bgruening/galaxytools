"""
Checks for accepting $schemas directive

run individually as py.test -k tests/test_schemas_directive.py
"""

import os
from typing import Any, Optional, Union

from schema_salad.avro.schema import Names, SchemaParseException
from schema_salad.ref_resolver import Loader
from schema_salad.schema import load_and_validate, load_schema

from .util import get_data

test_dir_name = "tests/"


class TestSchemasDirective:
    """Ensure codegen-produced parsers accept $schemas directives"""

    document_loader: Loader
    avsc_names: Union[Names, SchemaParseException, None] = None
    schema_metadata: Optional[dict[str, Any]] = None
    metaschema_loader: Optional[Loader] = None

    @classmethod
    def setup_class(cls) -> None:
        path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
        (
            cls.document_loader,
            cls.avsc_names,
            schema_metadata,
            metaschema_loader,
        ) = load_schema(path)

    def load_cwl(self, src: str) -> tuple[Any, dict[str, Any]]:
        path = get_data(test_dir_name + src)
        assert isinstance(self.avsc_names, Names)
        res = load_and_validate(self.document_loader, self.avsc_names, path, True)
        return res

    def test_dollarsign_schema(self) -> None:
        """EDAM.owl as a schema"""
        res = self.load_cwl(src="formattest2.cwl")

        # EDAM.owl resides in this directory
        assert os.path.split(str(res[0]["$schemas"][0]))[1] == "EDAM.owl"
