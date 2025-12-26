# SPDX-License-Identifier: Apache-2.0
"""
CWL file formats utilities.

For more information, please visit https://www.commonwl.org/user_guide/16-file-formats/
"""

from typing import Optional, Union

from rdflib import OWL, RDFS, Graph, URIRef
from schema_salad.exceptions import ValidationException
from schema_salad.utils import aslist, json_dumps

from cwl_utils.types import CWLObjectType


def formatSubclassOf(
    fmt: str, cls: str, ontology: Optional[Graph], visited: set[str]
) -> bool:
    """Determine if `fmt` is a subclass of `cls`."""
    if URIRef(fmt) == URIRef(cls):
        return True

    if ontology is None:
        return False

    if fmt in visited:
        return False

    visited.add(fmt)

    uriRefFmt = URIRef(fmt)

    for _s, _p, o in ontology.triples((uriRefFmt, RDFS.subClassOf, None)):
        # Find parent classes of `fmt` and search upward
        if formatSubclassOf(o, cls, ontology, visited):
            return True

    for _s, _p, o in ontology.triples((uriRefFmt, OWL.equivalentClass, None)):
        # Find equivalent classes of `fmt` and search horizontally
        if formatSubclassOf(o, cls, ontology, visited):
            return True

    for s, _p, _o in ontology.triples((None, OWL.equivalentClass, uriRefFmt)):
        # Find equivalent classes of `fmt` and search horizontally
        if formatSubclassOf(s, cls, ontology, visited):
            return True

    return False


def check_format(
    actual_file: Union[CWLObjectType, list[CWLObjectType]],
    input_formats: Union[list[str], str],
    ontology: Optional[Graph],
) -> None:
    """Confirm that the format present is valid for the allowed formats."""
    for afile in aslist(actual_file):
        if not afile:
            continue
        if "format" not in afile:
            raise ValidationException(
                f"File has no 'format' defined: {json_dumps(afile, indent=4)}"
            )
        for inpf in aslist(input_formats):
            if afile["format"] == inpf or formatSubclassOf(
                afile["format"], inpf, ontology, set()
            ):
                return
        raise ValidationException(
            f"File has an incompatible format: {json_dumps(afile, indent=4)}"
        )
