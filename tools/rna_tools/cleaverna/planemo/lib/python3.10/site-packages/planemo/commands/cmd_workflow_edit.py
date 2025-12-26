"""Module describing the planemo ``workflow_edit`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.galaxy.serve import sleep_for_serve
from planemo.runnable_resolve import for_runnable_identifier


@click.command("workflow_edit")
@options.required_workflow_arg()
@options.galaxy_serve_options()
@command_function
def cli(ctx, workflow_identifier, output=None, force=False, **kwds):
    """Open a synchronized Galaxy workflow editor."""
    assert is_galaxy_engine(**kwds)
    runnable = for_runnable_identifier(ctx, workflow_identifier, kwds)

    kwds["workflows_from_path"] = True

    with engine_context(ctx, **kwds) as galaxy_engine:
        with galaxy_engine.ensure_runnables_served([runnable]) as config:
            workflow_id = config.workflow_id_for_runnable(runnable)
            url = f"{config.galaxy_url}/workflows/edit?id={workflow_id}"
            click.launch(url)
            if kwds["engine"] != "external_galaxy":
                sleep_for_serve()
