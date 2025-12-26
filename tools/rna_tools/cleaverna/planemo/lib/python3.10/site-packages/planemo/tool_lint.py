from os.path import basename
from typing import (
    Any,
    Dict,
    TYPE_CHECKING,
)

from galaxy.tool_util.lint import lint_tool_source

import planemo.linters.biocontainer_registered
import planemo.linters.conda_requirements
import planemo.linters.doi
import planemo.linters.urls
import planemo.linters.xsd
from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK,
)
from planemo.io import (
    coalesce_return_codes,
    error,
    info,
)
from planemo.lint import build_lint_args
from planemo.tools import (
    is_tool_load_error,
    yield_tool_sources_on_paths,
)

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext

LINTING_TOOL_MESSAGE = "Linting tool %s"


def build_tool_lint_args(ctx: "PlanemoCliContext", **kwds) -> Dict[str, Any]:
    lint_args = build_lint_args(ctx, **kwds)
    extra_modules = _lint_extra_modules(**kwds)
    lint_args["extra_modules"] = extra_modules
    return lint_args


def lint_tools_on_path(ctx, paths, lint_args, **kwds):
    assert_tools = kwds.get("assert_tools", True)
    recursive = kwds.get("recursive", False)
    exit_codes = []
    for tool_path, tool_xml in yield_tool_sources_on_paths(ctx, paths, recursive):
        if handle_tool_load_error(tool_path, tool_xml):
            exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
            continue
        info("Linting tool %s" % tool_path)
        if not lint_tool_source(tool_xml, name=basename(tool_path), **lint_args):
            error("Failed linting")
            exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
        else:
            exit_codes.append(EXIT_CODE_OK)
    return coalesce_return_codes(exit_codes, assert_at_least_one=assert_tools)


def _lint_extra_modules(**kwds):
    linters = []

    if kwds.get("doi", False):
        linters.append(planemo.linters.doi)

    if kwds.get("urls", False):
        linters.append(planemo.linters.urls)

    if kwds.get("conda_requirements", False):
        linters.append(planemo.linters.conda_requirements)

    if kwds.get("biocontainer", False):
        linters.append(planemo.linters.biocontainer_registered)

    return linters


def handle_tool_load_error(tool_path, tool_xml):
    """Return True if tool_xml is tool load error (invalid XML), and
    print a helpful error message.
    """
    is_error = False
    if is_tool_load_error(tool_xml):
        info("Could not lint %s due to malformed xml." % tool_path)
        is_error = True
    return is_error
