"""Module describing the planemo ``dockstore_init`` command."""

import os

import click

from planemo import options
from planemo.cli import command_function
from planemo.io import launch_if_open_flagged
from planemo.workflow_lint import (
    DOCKSTORE_REGISTRY_CONF,
    generate_dockstore_yaml,
)


@click.command("dockstore_init")
@options.optional_project_arg()
@options.publish_dockstore_option()
@options.open_file_option()
@command_function
def cli(ctx, path=".", **kwds):
    """Initialize a .dockstore.yml configuration file for workflows in directory.

    Walk supplied directory and find all Galaxy workflows and test configurations
    and create a ``.dockstore.yml`` with references to these files. Be sure to push
    this file to Github before registering your workflow repository with Dockstore.

    Visit Dockstore at https://dockstore.org/. Find more about registering workflows
    with Dockstore at
    https://docs.dockstore.org/en/develop/getting-started/dockstore-workflows.html.
    """
    # TODO: implement -f semantics
    dockstore_path = os.path.join(path, DOCKSTORE_REGISTRY_CONF)
    contents = generate_dockstore_yaml(path, kwds["publish"])
    with open(dockstore_path, "w") as f:
        f.write(contents)
    launch_if_open_flagged(dockstore_path, **kwds)
