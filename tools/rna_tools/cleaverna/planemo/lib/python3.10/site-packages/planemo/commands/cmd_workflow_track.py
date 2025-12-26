"""Module describing the planemo ``workflow_track`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine.factory import engine_context
from planemo.galaxy.activity import wait_for_invocation_and_jobs


@click.command("workflow_track")
@options.invocation_target_options()
@options.fail_fast_option()
@command_function
def cli(ctx, invocation_id, **kwds):
    """Follow the progress of a workflow invocation."""
    with engine_context(ctx, engine="external_galaxy", **kwds) as engine, engine.ensure_runnables_served([]) as config:
        user_gi = config.user_gi
        wait_for_invocation_and_jobs(
            ctx,
            invocation_id,
            history_id=None,
            user_gi=user_gi,
            polling_backoff=5,
            fail_fast=kwds.get("fail_fast", False),
        )

    ctx.exit(0)
