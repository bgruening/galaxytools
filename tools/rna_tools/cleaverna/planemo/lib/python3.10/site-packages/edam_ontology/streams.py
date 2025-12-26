"""Return text streams to EDAM data.

Keep parsing and such to minimum to maintain maximum backward compatibility.
"""
import importlib.resources
from typing import IO


def tabular_stream() -> IO[str]:
    """Yield EDAM data in TSV format as a Python UTF-8 encoded text stream."""
    ref = importlib.resources.files("edam_ontology").joinpath("EDAM.tsv")
    return ref.open()
