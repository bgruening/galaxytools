"""Module describing the planemo ``conda_build`` command."""

from typing import (
    Tuple,
    TYPE_CHECKING,
)

import click

from planemo import options
from planemo.cli import command_function
from planemo.conda import build_conda_context
from planemo.io import error

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext


@click.command("conda_build")
@options.conda_target_options(include_local=False)  # No reason to expose local, we have to use it.
@options.recipe_arg(multiple=True)
@command_function
def cli(ctx: "PlanemoCliContext", paths: Tuple[str], **kwds) -> None:
    """Perform conda build with Planemo's conda."""
    # Force conda_use_local for building...
    kwds["conda_use_local"] = True
    conda_context = build_conda_context(ctx, handle_auto_init=True, **kwds)
    exit_code = conda_context.exec_command("build", paths)
    if exit_code:
        error(f"Failed to build [{' '.join(paths)}] with conda.")
    ctx.exit(exit_code)
