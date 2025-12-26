"""Module describing the planemo ``shed_create`` command."""

import sys

import click

from planemo import (
    options,
    shed,
)
from planemo.cli import command_function
from planemo.io import info


@click.command("shed_create")
@options.shed_publish_options()
@options.shed_message_option()
@options.shed_skip_upload()
@command_function
def cli(ctx, paths, **kwds):
    """Create a repository in a Galaxy Tool Shed.

    This will read the settings from the ``.shed.yml`` file.
    """
    shed_context = shed.get_shed_context(ctx, **kwds)

    def create(realized_repository):
        repo_id = realized_repository.find_repository_id(ctx, shed_context)
        if repo_id is None:
            if realized_repository.create(ctx, shed_context):
                info("Repository created")
                if not kwds["skip_upload"]:
                    return shed.upload_repository(ctx, realized_repository, **kwds)
                else:
                    return 0
            else:
                return 2
        else:
            return 1

    exit_code = shed.for_each_repository(ctx, create, paths, **kwds)
    sys.exit(exit_code)
