"""Module describing the planemo ``conda_init`` command."""

import click
from galaxy.tool_util.deps import conda_util

from planemo import options
from planemo.cli import command_function
from planemo.conda import build_conda_context
from planemo.exit_codes import EXIT_CODE_ALREADY_EXISTS
from planemo.io import (
    info,
    warn,
)

MESSAGE_ERROR_ALREADY_EXISTS = "conda_init failed - Conda appears to already be installed at '%s'"
MESSAGE_ERROR_FAILED = "conda_init failed - failed to install to '%s'"
MESSAGE_INSTALL_OKAY = "Conda installation succeeded - Conda is available at '%s'"


@click.command("conda_init")
@options.conda_target_options(include_local=False)  # Always use local during init.
@command_function
def cli(ctx, **kwds):
    """Download and install conda.

    This will download conda for managing dependencies for your platform
    using the appropriate Miniconda installer.

    By running this command, you are agreeing to the terms of the conda
    license a 3-clause BSD 3 license. Please review full license at
    http://docs.continuum.io/anaconda/eula.

    Planemo will print a warning and terminate with an exit code of 7
    if Conda is already installed.
    """
    conda_context = build_conda_context(ctx, **kwds)
    if conda_context.is_conda_installed():
        warn(MESSAGE_ERROR_ALREADY_EXISTS % conda_context.conda_exec)
        exit = EXIT_CODE_ALREADY_EXISTS
    else:
        exit = conda_util.install_conda(conda_context=conda_context, force_conda_build=True)
        if exit:
            warn(MESSAGE_ERROR_FAILED % conda_context.conda_exec)
        else:
            info(MESSAGE_INSTALL_OKAY % conda_context.conda_exec)

    ctx.exit(exit)
