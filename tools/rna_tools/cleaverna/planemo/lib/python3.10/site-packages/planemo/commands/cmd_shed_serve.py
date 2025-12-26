"""Module describing the planemo ``shed_serve`` command."""

import click

from planemo import (
    io,
    options,
    shed,
)
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.galaxy.serve import sleep_for_serve
from planemo.runnable_resolve import install_args_list_to_runnables


@click.command("shed_serve")
@options.shed_read_options()
@options.galaxy_serve_options()
@click.option(
    "--skip_dependencies", is_flag=True, help="Do not install shed dependencies as part of repository installation."
)
@command_function
def cli(ctx, paths, **kwds):
    """Launch Galaxy with Tool Shed dependencies.

    This command will start a Galaxy instance configured to target the
    specified shed, find published artifacts (tools and dependencies)
    corresponding to command-line arguments and ``.shed.yml`` file(s),
    install these artifacts, and serve a Galaxy instances that can be
    logged into and explored interactively.
    """
    install_args_list = kwds["install_args_list"] = shed.install_arg_lists(ctx, paths, **kwds)
    runnables = install_args_list_to_runnables(ctx, install_args_list, kwds)
    with engine_context(ctx, **kwds) as engine:
        with engine.ensure_runnables_served(runnables) as config:
            io.info(f"Galaxy running with tools installed at {config.galaxy_url}")
            sleep_for_serve()
