"""Module describing the planemo ``tool_init`` command."""

from typing import (
    Any,
    Dict,
)

import click

from planemo import (
    io,
    options,
    tool_builder,
)
from planemo.cli import command_function
from planemo.options import tool_init_autopygen_option


@click.command("tool_init")
@options.tool_init_id_option()
@options.force_option(what="tool")
@options.tool_init_tool_option()
@options.tool_init_name_option()
@options.tool_init_version_option()
@options.tool_init_description_option()
@options.tool_init_command_option()
@options.tool_init_example_command_option()
@options.tool_init_example_input_option()
@options.tool_init_example_output_option()
@options.tool_init_named_output_option()
@options.tool_init_input_option()
@options.tool_init_output_option()
@options.tool_init_help_text_option()
@options.tool_init_help_from_command_option()
@options.tool_init_doi_option()
@options.tool_init_cite_url_option()
@options.tool_init_test_case_option()
@options.tool_init_macros_option()
@options.tool_init_version_command_option()
@options.tool_init_requirement_option()
@options.tool_init_container_option()
@options.build_cwl_option()
@tool_init_autopygen_option()
@command_function
def cli(ctx, **kwds):
    """Generate tool outline from given arguments."""
    invalid = _validate_kwds(kwds)
    if invalid:
        ctx.exit(invalid)
    tool_description = tool_builder.build(**kwds)
    tool_builder.write_tool_description(ctx, tool_description, **kwds)


def _validate_kwds(kwds: Dict[str, Any]) -> int:
    def not_exclusive(x, y):
        if kwds.get(x) and kwds.get(y):
            io.error(f"Can only specify one of --{x} and --{y}")
            return True

    def not_specifing_dependent_option(x, y):
        if kwds.get(x) and not kwds.get(y):
            template = "Can only use the --%s option if also specifying --%s"
            message = template % (x, y)
            io.error(message)
            return True

    if not_exclusive("help_text", "help_from_command"):
        return 1
    if not_exclusive("command", "example_command"):
        return 1
    if not_specifing_dependent_option("example_input", "example_command"):
        return 1
    if not_specifing_dependent_option("example_output", "example_command"):
        return 1
    if not_specifing_dependent_option("test_case", "example_command"):
        return 1
    return 0
