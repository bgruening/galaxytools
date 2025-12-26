"""
Checks loading of some real world tools and workflows found in the wild (e.g. dockstore)

run individually as py.test -k tests/test_real_cwl.py
"""

from typing import Any, Optional, Union

import pytest

from schema_salad.avro.schema import Names, SchemaParseException
from schema_salad.exceptions import ValidationException
from schema_salad.ref_resolver import Loader
from schema_salad.schema import load_and_validate, load_schema

from .util import get_data

test_dir_name = "tests/test_real_cwl/"


class TestRealWorldCWL:
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

    def load_cwl(self, src: str) -> None:
        path = get_data(test_dir_name + src)
        assert isinstance(self.avsc_names, Names)
        with pytest.raises(ValidationException):
            try:
                load_and_validate(
                    self.document_loader,
                    self.avsc_names,
                    path,
                    True,
                )
            except ValidationException as e:
                # msgs = to_one_line_messages(str(e)).splitlines()
                print("\n", e)
                raise

    def test_topmed_single_doc(self) -> None:
        """TOPMed Variant Calling Pipeline CWL1"""
        self.load_cwl(src="topmed/topmed_variant_calling_pipeline.cwl")

    def test_h3agatk_WES(self) -> None:
        """H3ABioNet GATK Germline Workflow"""
        self.load_cwl(src="h3agatk/GATK-complete-WES-Workflow-h3abionet.cwl")

    def test_h3agatk_SNP(self) -> None:
        """H3ABioNet SNPs Workflow"""
        self.load_cwl(src="h3agatk/GATK-Sub-Workflow-h3abionet-snp.cwl")

    def test_icgc_pancan(self) -> None:
        """ICGC PanCan"""
        self.load_cwl(src="ICGC-TCGA-PanCancer/preprocess_vcf.cwl")
