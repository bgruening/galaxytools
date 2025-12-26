"""Module contains factory method for building class:`Engine` objects."""

import contextlib
from typing import Generator

from planemo.engine.interface import BaseEngine
from .cwltool import CwlToolEngine
from .galaxy import (
    DockerizedManagedGalaxyEngine,
    ExternalGalaxyEngine,
    LocalManagedGalaxyEngine,
    LocalManagedGalaxyEngineWithSingularityDB,
)
from .toil import ToilEngine

UNKNOWN_ENGINE_TYPE_MESSAGE = "Unknown engine type specified [%s]."


def is_galaxy_engine(**kwds):
    """Return True iff the engine configured is :class:`GalaxyEngine`."""
    engine_type_str = kwds.get("engine", "galaxy")
    return engine_type_str in ["galaxy", "docker_galaxy", "external_galaxy"]


def build_engine(ctx, **kwds):
    """Build an engine from the supplied planemo configuration."""
    engine_type_str = kwds.get("engine", "galaxy")
    if engine_type_str == "galaxy":
        if "database_type" in kwds and kwds["database_type"] == "postgres_singularity":
            engine_type = LocalManagedGalaxyEngineWithSingularityDB
        else:
            engine_type = LocalManagedGalaxyEngine
    elif engine_type_str == "docker_galaxy":
        engine_type = DockerizedManagedGalaxyEngine
    elif engine_type_str == "external_galaxy":
        engine_type = ExternalGalaxyEngine
    elif engine_type_str == "cwltool":
        engine_type = CwlToolEngine
    elif engine_type_str == "toil":
        engine_type = ToilEngine
    else:
        raise Exception(UNKNOWN_ENGINE_TYPE_MESSAGE % engine_type_str)

    return engine_type(ctx, **kwds)


@contextlib.contextmanager
def engine_context(ctx, **kwds) -> Generator[BaseEngine, None, None]:
    """A :func:`contextlib.contextmanager` engine builder for use with ``with`` statements.

    https://docs.python.org/2/library/contextlib.html
    """
    engine = None
    try:
        engine = build_engine(ctx, **kwds)
        yield engine
    finally:
        if engine is not None:
            engine.cleanup()


__all__ = (
    "is_galaxy_engine",
    "build_engine",
    "engine_context",
)
