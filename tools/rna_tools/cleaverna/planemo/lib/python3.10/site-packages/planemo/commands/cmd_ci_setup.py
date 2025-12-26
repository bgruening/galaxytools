"""Module describing the planemo ``ci_setup`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.serve import serve_daemon


@click.command("ci_setup")
@options.galaxy_target_options()
@command_function
def cli(ctx, **kwds):
    """
    Launch Galaxy instance, then terminate instance.

    Useful for populating a CI cache.
    """
    kwds["galaxy_skip_client_build"] = True
    kwds["no_dependency_resolution"] = True
    with serve_daemon(ctx, **kwds):
        return
