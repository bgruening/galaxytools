"""Module describing the planemo ``project_init`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.io import (
    info,
    launch_if_open_flagged,
    warn,
)

SUCCESS_MESSAGE = "Wrote configuration template to %s for your profile, please open with editor to customize if needed."


@click.command("profile_job_config_init")
@options.profile_name_argument()
@options.open_file_option()
@options.job_config_init_options()
@command_function
def cli(ctx, profile_name, **kwds):
    """Initialize a Galaxy job configuration for specified profile."""
    if not profiles.profile_exists(ctx, profile_name, **kwds):
        warn(f"Profile '{profile_name}' does not exist, exiting.")
        pass

    config_path = profiles.initialize_job_config(ctx, profile_name, **kwds)
    info(SUCCESS_MESSAGE % config_path)
    launch_if_open_flagged(config_path, **kwds)
