"""Module describing the planemo ``profile_create`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles


@click.command("profile_create")
@options.profile_name_argument()
@options.profile_database_options()
@options.serve_engine_option()
@options.docker_config_options()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@options.galaxy_admin_key_option()
@command_function
def cli(ctx, profile_name, **kwds):
    """Create a profile."""
    profiles.create_profile(ctx, profile_name, **kwds)
    print("Profile [%s] created." % profile_name)
