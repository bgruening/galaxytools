"""Module describing the planemo ``syntax`` command."""

import click

from planemo.cli import command_function

SYNTAX_URL = "https://docs.galaxyproject.org/en/latest/dev/schema.html"


@click.command("syntax")
@command_function
def cli(ctx, **kwds):
    """Open tool config syntax page in web browser."""
    click.launch(SYNTAX_URL)
