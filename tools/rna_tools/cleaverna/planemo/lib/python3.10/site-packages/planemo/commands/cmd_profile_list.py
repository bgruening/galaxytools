"""Module describing the planemo ``profile_list`` command."""

import click

from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.io import info


@click.command("profile_list")
@command_function
def cli(ctx, **kwds):
    """List configured profile names."""
    info("Looking for profiles...")
    profile_names = profiles.list_profiles(ctx, **kwds)
    for profile in profile_names:
        print(profile)
    info(f"{len(profile_names)} configured profiles are available.")
