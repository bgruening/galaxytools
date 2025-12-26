"""Module describing the planemo ``training_generate_from_wf`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.training import Training


@click.command("training_generate_from_wf")
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.training_generate_tuto_from_wf_options()
@options.galaxy_serve_options()
@command_function
def cli(ctx, uris, **kwds):
    """Create tutorial skeleton from workflow."""
    kwds["no_dependency_resolution"] = True
    # Ugh rather than override this - it just needs to not be in the serve options
    # for this command.
    kwds["skip_client_build"] = True
    training = Training(kwds)
    training.generate_tuto_from_wf(ctx)
