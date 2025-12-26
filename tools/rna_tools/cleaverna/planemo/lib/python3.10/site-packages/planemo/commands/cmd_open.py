"""Module describing the planemo ``open`` command."""

import click

from planemo.cli import command_function


@click.command("open")
@click.argument(
    "path",
    metavar="PATH",
    type=click.Path(
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        resolve_path=True,
    ),
    default="tool_test_output.html",
)
@command_function
def cli(ctx, path, **kwds):
    """Open latest Planemo test results in a web browser."""
    click.launch(path)
