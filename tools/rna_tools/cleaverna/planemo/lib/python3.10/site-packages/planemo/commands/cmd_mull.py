"""Module describing the planemo ``mull`` command."""

import click
from galaxy.tool_util.deps.mulled.mulled_build import mull_targets

from planemo import options
from planemo.cli import command_function
from planemo.mulled import (
    build_mull_target_kwds,
    collect_mulled_target_lists,
)


@click.command("mull")
@options.optional_tools_arg(multiple=True)
@options.recursive_option()
@options.mulled_options()
@options.conda_ensure_channels_option()
@command_function
def cli(ctx, paths, **kwds):
    """Build containers for specified tools.

    Supplied tools will be inspected for referenced requirement packages. For
    each combination of requirements a "mulled" container will be built. Galaxy
    can automatically discover this container and subsequently use it to run
    or test the tool.

    For this to work, the tool's requirements will need to be present in a known
    Conda channel such as bioconda (https://github.com/bioconda/bioconda-recipes).
    This can be verified by running ``planemo lint --conda_requirements`` on the
    target tool(s).
    """
    for mulled_targets in collect_mulled_target_lists(ctx, paths, recursive=kwds["recursive"]):
        mull_target_kwds = build_mull_target_kwds(ctx, **kwds)
        command = kwds["mulled_command"]
        mull_targets(mulled_targets, command=command, **mull_target_kwds)
