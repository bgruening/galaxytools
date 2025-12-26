"""Module describing :class:`planemo.interface.Engine` abstraction and implementations.

The main entrypoint is the ``contextmanager`` :func:`engine_context`.
"""

from .factory import (
    engine_context,
    is_galaxy_engine,
)

__all__ = (
    "engine_context",
    "is_galaxy_engine",
)
