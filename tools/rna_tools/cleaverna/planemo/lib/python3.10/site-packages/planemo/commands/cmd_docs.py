"""Module describing the planemo ``docs`` command."""

import click

from planemo.cli import command_function

DOCS_URL = "http://planemo.readthedocs.io/en/latest/"


@click.command("docs")
@command_function
def cli(ctx, **kwds):
    """Open Planemo documentation in web browser."""
    click.launch(DOCS_URL)
