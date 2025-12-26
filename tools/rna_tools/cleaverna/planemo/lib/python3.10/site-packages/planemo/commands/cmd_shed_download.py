"""Module describing the planemo ``shed_download`` command."""

import sys

import click

from planemo import (
    options,
    shed,
)
from planemo.cli import command_function

target_path = click.Path(
    file_okay=True,
    writable=True,
    resolve_path=True,
)


@click.command("shed_download")
@options.shed_read_options()
@click.option(
    "--destination",
    default="shed_download.tar.gz",
    type=target_path,
    help="Destination pattern of tarball(s) to download - if this doesn't "
    "end in 'gz' it will be treated as a directory to extract tool "
    "contents into (defaults to shed_download.tar.gz). If multiple "
    "repositories are discovered in a .shed.yml file these will be "
    "created as shed_download_<name>.tar.gz by default for instance, "
    "simpler repositories will just be downloaded to the specified file.",
)
@command_function
def cli(ctx, paths, **kwds):
    """Download tool from Tool Shed into directory.

    Download a tool repository as a tarball from the tool shed and extract
    to the specified directory.
    """
    shed_context = shed.get_shed_context(ctx, read_only=True, **kwds)

    def download(realized_repository):
        return shed.download_tarball(ctx, shed_context, realized_repository, **kwds)

    exit_code = shed.for_each_repository(ctx, download, paths, **kwds)
    sys.exit(exit_code)
