"""Module describing the planemo ``mulled_init`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.mulled import build_involucro_context


@click.command("mulled_init")
@options.mulled_options()
@command_function
def cli(ctx, **kwds):
    """Download and install involucro for mull command.

    This will happen automatically when using the mull command, but this can
    be pre-installed in an environment using this command.
    """
    build_involucro_context(ctx, **kwds)
