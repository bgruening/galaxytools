"""Module describing the planemo ``shed_upload`` command."""

import sys

import click

from planemo import (
    options,
    shed,
)
from planemo.cli import command_function

tar_path = click.Path(
    exists=True,
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
)


@click.command("shed_upload")
@options.shed_publish_options()
@options.shed_upload_options()
@click.option(
    "--tar_only",
    is_flag=True,
    help="Produce tar file for upload but do not publish to a tool shed.",
)
@click.option(
    "--tar",
    help="Specify a pre-existing tar file instead of automatically building one as part of this command.",
    type=tar_path,
    default=None,
)
@command_function
def cli(ctx, paths, **kwds):
    """Low-level command to upload tarballs.

    Generally, ``shed_update`` should be used instead since it also updates
    both tool shed contents (via tar ball generation and upload) as well as
    metadata (to handle metadata changes in ``.shed.yml`` files).

    \b
        % planemo shed_upload --tar_only  ~/
        % tar -tzf shed_upload.tar.gz
        test-data/blastdb.loc
        ...
        tools/ncbi_blast_plus/tool_dependencies.xml
        % tar -tzf shed_upload.tar.gz | wc -l
        117

    """

    def upload(realized_repository):
        return shed.upload_repository(ctx, realized_repository, **kwds)

    exit_code = shed.for_each_repository(ctx, upload, paths, **kwds)
    sys.exit(exit_code)
