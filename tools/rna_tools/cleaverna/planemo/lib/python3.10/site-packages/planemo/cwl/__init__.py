"""Entry point for modules describing abstractions for dealing with CWL artifacts."""

from .run import run_cwltool
from .toil import run_toil

__all__ = (
    "run_cwltool",
    "run_toil",
)
