"""Module describing the planemo ``shed_build`` command."""

import shutil
import sys

import click

from planemo import (
    options,
    shed,
)
from planemo.cli import command_function


@click.command("shed_build")
@options.optional_tools_arg(multiple=False)
@command_function
def cli(ctx, path, **kwds):
    """Create a Galaxy tool tarball.

    This will use the .shed.yml file to prepare a tarball
    (which you could upload to the Tool Shed manually).
    """

    def build(realized_repository):
        tarpath = shed.build_tarball(realized_repository.path)
        outpath = realized_repository.real_path + ".tar.gz"
        shutil.move(tarpath, outpath)
        print("Created: %s" % (outpath))
        return 0

    exit_code = shed.for_each_repository(ctx, build, [path], **kwds)
    sys.exit(exit_code)
