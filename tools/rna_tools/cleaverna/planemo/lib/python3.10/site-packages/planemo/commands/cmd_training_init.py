"""Module describing the planemo ``training_init`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.training import Training


@click.command("training_init")
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.training_init_options()
# @options.force_option()
@options.galaxy_serve_options()
@command_function
def cli(ctx, uris, **kwds):
    """Build training template from workflow."""
    kwds["no_dependency_resolution"] = True
    training = Training(kwds)
    training.init_training(ctx)
