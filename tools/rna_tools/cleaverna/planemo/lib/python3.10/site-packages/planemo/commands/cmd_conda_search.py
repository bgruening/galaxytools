"""Module describing the planemo ``conda_search`` command."""

import click
import packaging.version

from planemo import options
from planemo.cli import command_function
from planemo.conda import build_conda_context

VERSION_4_DOT_4 = packaging.version.Version("4.4")


@click.command("conda_search")
@options.conda_target_options(include_local=False)
@click.argument(
    "term",
    metavar="TERM",
    type=str,
    nargs=1,
)
@command_function
def cli(ctx, term, **kwds):
    """Perform conda search with Planemo's conda.

    Implicitly adds channels Planemo is configured with.
    """
    conda_context = build_conda_context(ctx, handle_auto_init=True, **kwds)
    # Handle CLI interface change for conda search in 4.5.
    #  xref: https://github.com/conda/conda/pull/5597/files
    if conda_context.conda_version >= VERSION_4_DOT_4:
        term = "*%s*" % term
    args = conda_context._override_channels_args + [term]
    exit_code = conda_context.exec_command("search", args)
    ctx.exit(exit_code)
