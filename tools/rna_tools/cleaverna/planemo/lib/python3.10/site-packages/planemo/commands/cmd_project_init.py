"""Module describing the planemo ``project_init`` command."""

import os
import shutil
import tempfile

import click

from planemo import options
from planemo.cli import command_function
from planemo.io import (
    untar_to,
    warn,
)

SOURCE_HOST = "https://codeload.github.com"
DOWNLOAD_URL = "%s/galaxyproject/planemo/tar.gz/master" % SOURCE_HOST


@click.command("project_init")
@options.optional_project_arg(exists=None)
@click.option("--template", default=None)
@command_function
def cli(ctx, path, template=None, **kwds):
    """(Experimental) Initialize a new tool project.

    This is only a proof-of-concept demo right now.
    """
    if template is None:
        warn("Creating empty project, this function doesn't do much yet.")
    if not os.path.exists(path):
        os.makedirs(path)
    if template is None:
        return

    tempdir = tempfile.mkdtemp()
    tar_args = ["-zxf", "-", "--strip-components=2"]
    try:
        untar_to(DOWNLOAD_URL, tar_args=tar_args, dest_dir=tempdir)
        template_dir = os.path.join(tempdir, template)
        for entry in os.listdir(template_dir):
            shutil.move(os.path.join(template_dir, entry), path)
    finally:
        shutil.rmtree(tempdir)
