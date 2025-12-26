"""Module describing the planemo ``rerun`` command."""

from typing import (
    Tuple,
    TYPE_CHECKING,
)

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.engine.galaxy import ExternalGalaxyEngine
from planemo.io import (
    error,
    info,
)
from planemo.runnable import Rerunnable

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext


@click.command("rerun")
@options.profile_option()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@click.option(
    "--invocation",
    "rerunnable_type",
    flag_value="invocation",
    help=("Rerun failed jobs associated by one or more invocation IDs."),
)
@click.option(
    "--history",
    "rerunnable_type",
    flag_value="history",
    help=("Rerun failed jobs associated by one or more history IDs."),
)
@click.option(
    "--job", "rerunnable_type", flag_value="job", help=("Rerun failed jobs specified by one or more job IDs.")
)
@click.argument(
    "rerunnable_ids",
    metavar="RERUNNABLE_IDS",
    type=str,
    nargs=-1,
)
@command_function
def cli(ctx: "PlanemoCliContext", rerunnable_ids: Tuple[str], **kwds) -> None:
    """Planemo command for rerunning and remapping failed jobs on an external Galaxy server.
    Supply a list of history, invocation or job IDs, identifying the ID type using the
    --invocation, --history or --job flag, and all associated failed jobs will be rerun.

    Please note: attempting to rerun non-remappable jobs will result in an exit code of 1. As
    jobs cannot be remapped more than once, running this command two or more times with the same
    history or job IDs will therefore return an exit code of 1. If avoiding this is important,
    you should specify the invocation ID instead if possible.

    \b
        % planemo rerun --invocation / --history / --job RERUNNABLE_IDS
    """
    # Possible TODO: allow collection IDs to be specified as well
    if not kwds.get("rerunnable_type"):
        error("Please specify the type (invocation, history or job) of the IDs which should be rerun.")
        ctx.exit(1)
    kwds["engine"] = "external_galaxy"
    rerun_successful = True
    with engine_context(ctx, **kwds) as engine:
        assert isinstance(engine, ExternalGalaxyEngine)
        for rerunnable_id in rerunnable_ids:
            rerunnable = Rerunnable(rerunnable_id, kwds["rerunnable_type"], kwds["galaxy_url"])
            rerun_result = engine.rerun(ctx, rerunnable, **kwds)
            if not rerun_result.was_successful:
                rerun_successful = False

    if rerun_successful:
        info("All requested jobs were rerun successfully.")
        ctx.exit(0)
    else:
        error("Some of the requested jobs could not be rerun.")
        ctx.exit(1)
