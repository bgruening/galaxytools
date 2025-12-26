"""Module describing the planemo ``profile_delete`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles


@click.command("profile_delete")
@options.profile_name_argument()
@options.profile_database_options()
@options.docker_config_options()
@command_function
def cli(ctx, profile_name, **kwds):
    """Delete a profile."""
    profiles.delete_profile(ctx, profile_name, **kwds)
    print("Profile deleted.")
