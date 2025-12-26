"""Planemo specific utilities for dealing with mulled containers.

The extend galaxy-tool-util's features with planemo specific idioms.
"""

import os
from typing import (
    Iterable,
    List,
)

from galaxy.tool_util.deps.mulled.mulled_build import (
    DEFAULT_CHANNELS,
    ensure_installed,
    InvolucroContext,
)
from galaxy.tool_util.deps.mulled.util import (
    build_target,
    CondaTarget,
)

from planemo.conda import collect_conda_target_lists
from planemo.io import shell


def conda_to_mulled_targets(conda_targets: Iterable[CondaTarget]) -> List[CondaTarget]:
    return list(map(lambda c: build_target(c.package, c.version), conda_targets))


def collect_mulled_target_lists(ctx, paths: Iterable[str], recursive: bool = False) -> List[List[CondaTarget]]:
    return list(map(conda_to_mulled_targets, collect_conda_target_lists(ctx, paths, recursive=recursive)))


def build_involucro_context(ctx, **kwds):
    """Build a galaxy-tool-util InvolucroContext tailored to planemo use.

    Using planemo's common command-line/global config options.
    """
    involucro_path_default = os.path.join(ctx.workspace, "involucro")
    involucro_path = kwds.get("involucro_path", involucro_path_default)
    use_planemo_shell = kwds.get("use_planemo_shell_exec", True)
    shell_exec = shell if use_planemo_shell else None
    involucro_context = InvolucroContext(involucro_bin=involucro_path, shell_exec=shell_exec)
    if not ensure_installed(involucro_context, True):
        raise Exception("Failed to install involucro for Planemo.")
    return involucro_context


def build_mull_target_kwds(ctx, **kwds):
    """Adapt Planemo's CLI and workspace configuration to galaxy-tool-util's mulled_build options."""
    involucro_context = build_involucro_context(ctx, **kwds)
    channels = kwds.get("conda_ensure_channels", ",".join(DEFAULT_CHANNELS))
    namespace = kwds.get("mulled_namespace", "biocontainers")
    target_kwds = {
        "involucro_context": involucro_context,
        "channels": channels.split(","),
        "namespace": namespace,
        "verbose": ctx.verbose,
    }

    conda_version = kwds.get("mulled_conda_version", None)
    if conda_version is not None:
        target_kwds["conda_version"] = conda_version
    return target_kwds


__all__ = (
    "build_involucro_context",
    "build_mull_target_kwds",
    "collect_mulled_target_lists",
    "conda_to_mulled_targets",
)
