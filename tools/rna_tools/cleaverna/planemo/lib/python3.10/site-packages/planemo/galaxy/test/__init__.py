"""Entry point and interface for ``planemo.galaxy.test`` package."""

from .actions import (
    handle_reports,
    handle_reports_and_summary,
)
from .structures import StructuredData

__all__ = (
    "handle_reports",
    "handle_reports_and_summary",
    "StructuredData",
)
