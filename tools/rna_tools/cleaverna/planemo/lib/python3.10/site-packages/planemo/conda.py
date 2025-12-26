"""Planemo specific utilities for dealing with conda.

The extend galaxy-tool-util's features with planemo specific idioms.
"""

import collections
import os
import threading
from copy import deepcopy
from typing import (
    Any,
    Dict,
    FrozenSet,
    Iterable,
    List,
    Optional,
    Set,
    Tuple,
    TYPE_CHECKING,
    Union,
)

from galaxy.tool_util.deps import conda_util
from galaxy.tool_util.deps.conda_util import (
    CondaContext,
    CondaTarget,
)
from galaxy.util import unicodify

from planemo.exit_codes import (
    EXIT_CODE_FAILED_DEPENDENCIES,
    ExitCodeException,
)
from planemo.io import (
    error,
    shell,
)
from planemo.tools import yield_tool_sources_on_paths

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

BEST_PRACTICE_CHANNELS = ["conda-forge", "bioconda"]


def build_conda_context(ctx: "PlanemoCliContext", **kwds) -> CondaContext:
    """Build a galaxy-tool-util CondaContext tailored to planemo use.

    Using planemo's common command-line/global config options.
    """
    condarc_override_default = os.path.join(ctx.workspace, "condarc")
    conda_prefix = kwds.get("conda_prefix", None)
    conda_exec = kwds.get("conda_exec", None)
    use_planemo_shell = kwds.get("use_planemo_shell_exec", True)
    ensure_channels = kwds.get("conda_ensure_channels", "")
    condarc_override = kwds.get("condarc", condarc_override_default)
    use_local = kwds.get("conda_use_local", False)
    shell_exec = shell if use_planemo_shell else None
    conda_context = CondaContext(
        conda_prefix=conda_prefix,
        ensure_channels=ensure_channels,
        condarc_override=condarc_override,
        use_local=use_local,
        shell_exec=shell_exec,
        conda_exec=conda_exec,
    )
    handle_auto_init = kwds.get("handle_auto_init", False)
    if handle_auto_init and not conda_context.is_installed():
        auto_init = kwds.get("conda_auto_init", True)
        failed = True
        if auto_init:
            if conda_context.can_install_conda():
                if conda_util.install_conda(conda_context) != 0:
                    error("Attempted to install conda and failed.")
                else:
                    failed = False
            else:
                error("Cannot install Conda - perhaps due to a failed installation or permission problems.")
        else:
            error("Conda not configured - run ``planemo conda_init`` or pass ``--conda_auto_init`` to continue.")

        if failed:
            raise ExitCodeException(EXIT_CODE_FAILED_DEPENDENCIES)
    if handle_auto_init:
        if conda_context.ensure_conda_build_installed_if_needed() != 0:
            error("Attempted to install conda-build and failed.")
    return conda_context


def collect_conda_targets(ctx, paths, recursive=False, found_tool_callback=None):
    """Load CondaTarget objects from supplied artifact sources.

    If a tool contains more than one requirement, the requirements will each
    appear once in the output.
    """
    conda_targets = set()
    real_paths = []
    for path in paths:
        if not os.path.exists(path):
            targets = target_str_to_targets(path)
            [conda_targets.add(_) for _ in targets]
        else:
            real_paths.append(path)

    for tool_path, tool_source in yield_tool_sources_on_paths(
        ctx, real_paths, recursive=recursive, exclude_deprecated=True
    ):
        if found_tool_callback:
            found_tool_callback(tool_path)
        for target in tool_source_conda_targets(tool_source):
            conda_targets.add(target)
    return conda_targets


# Copied and modified from mulled stuff - need to syncronize these concepts.
def target_str_to_targets(targets_raw: str) -> List[CondaTarget]:
    def parse_target(target_str: str) -> CondaTarget:
        if "=" in target_str:
            package_name, version = target_str.split("=", 1)
        else:
            package_name = target_str
            version = None
        target = CondaTarget(package_name, version)
        return target

    targets = [parse_target(_) for _ in targets_raw.split(",")]
    return targets


def collect_conda_target_lists(
    ctx: "PlanemoCliContext", paths: Iterable[str], recursive: bool = False, found_tool_callback=None
) -> List[FrozenSet[CondaTarget]]:
    """Load CondaTarget lists from supplied artifact sources.

    If a tool contains more than one requirement, the requirements will all
    appear together as one list element of the output list.
    """
    conda_target_lists, _ = collect_conda_target_lists_and_tool_paths(
        ctx, paths, recursive=recursive, found_tool_callback=found_tool_callback
    )
    return conda_target_lists


def collect_conda_target_lists_and_tool_paths(
    ctx: "PlanemoCliContext", paths: Iterable[str], recursive: bool = False, found_tool_callback=None
) -> Tuple[List[FrozenSet[CondaTarget]], List[List[str]]]:
    """Load CondaTarget lists from supplied artifact sources.

    If a tool contains more than one requirement, the requirements will all
    appear together as one list element of the output list.
    """
    conda_target_sets: Set[FrozenSet[CondaTarget]] = set()
    tool_paths = collections.defaultdict(list)
    for tool_path, tool_source in yield_tool_sources_on_paths(ctx, paths, recursive=recursive, yield_load_errors=False):
        try:
            if found_tool_callback:
                found_tool_callback(tool_path)
            targets = frozenset(tool_source_conda_targets(tool_source))
            conda_target_sets.add(targets)
            tool_paths[targets].append(tool_path)
        except Exception as e:
            ctx.log(f"Error while collecting list of conda targets for '{tool_path}': {unicodify(e)}")

    # Turn them into lists so the order matches before returning...
    conda_target_lists = list(conda_target_sets)
    conda_target_tool_paths = [tool_paths[c] for c in conda_target_lists]

    return conda_target_lists, conda_target_tool_paths


def tool_source_conda_targets(tool_source):
    """Load CondaTarget object from supplied abstract tool source."""
    requirements, *_ = tool_source.parse_requirements_and_containers()
    return conda_util.requirements_to_conda_targets(requirements)


best_practice_search_first = threading.local()


def best_practice_search(
    conda_target: CondaTarget, conda_context: Optional[CondaContext] = None, platform: Optional[str] = None
) -> Union[Tuple[None, None], Tuple[Dict[str, Any], bool]]:
    # Call it in offline mode after the first time.
    try:
        best_practice_search_first.previously_called
        # TODO: Undo this...
        offline = False
    except AttributeError:
        best_practice_search_first.previously_called = True
        offline = False

    if conda_context:
        if conda_context.ensure_channels != BEST_PRACTICE_CHANNELS:
            conda_context = deepcopy(conda_context)
            conda_context.ensure_channels = BEST_PRACTICE_CHANNELS
    else:
        conda_context = CondaContext(ensure_channels=BEST_PRACTICE_CHANNELS)
    return conda_util.best_search_result(
        conda_target,
        conda_context=conda_context,
        offline=offline,
        platform=platform,
    )


__all__ = (
    "BEST_PRACTICE_CHANNELS",
    "best_practice_search",
    "build_conda_context",
    "collect_conda_targets",
    "collect_conda_target_lists",
    "collect_conda_target_lists_and_tool_paths",
    "tool_source_conda_targets",
)
