"""Tool linting module that lints Galaxy tool against experimental XSD."""

import copy
import os
import tempfile

import galaxy.tool_util

import planemo.lint

TOOL_XSD = os.path.join(os.path.dirname(galaxy.tool_util.__file__), "xsd", "galaxy.xsd")


def lint_tool_xsd(tool_xml, lint_ctx):
    """Write a temp file out and lint it."""
    with tempfile.NamedTemporaryFile() as tf:
        _clean_root(tool_xml).write(tf.name)
        planemo.lint.lint_xsd(lint_ctx, TOOL_XSD, tf.name)


def _clean_root(tool_xml):
    """XSD assumes macros have been expanded, so remove them."""
    clean_tool_xml = copy.deepcopy(tool_xml)
    to_remove = []
    for macros_el in clean_tool_xml.getroot().findall("macros"):
        to_remove.append(macros_el)
    for macros_el in to_remove:
        clean_tool_xml.getroot().remove(macros_el)
    return clean_tool_xml
