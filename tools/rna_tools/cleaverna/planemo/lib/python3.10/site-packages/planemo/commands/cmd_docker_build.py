"""Module describing the planemo ``docker_build`` command."""

import click
from galaxy.tool_util.deps import dockerfiles

from planemo import options
from planemo.cli import command_function
from planemo.io import error


@click.command("docker_build")
@options.optional_tools_arg()
@click.option("--dockerfile", default=None)
@click.option("--docker_image_cache", default=None)
@options.docker_cmd_option()
@options.docker_sudo_option()
@options.docker_sudo_cmd_option()
@options.docker_host_option()
@command_function
def cli(ctx, path=".", dockerfile=None, **kwds):
    """Build (and optionally cache) Docker images.

    Loads the tool or tools referenced by ``TOOL_PATH`` (by default all tools
    in current directory), and ensures they all reference the same Docker image
    and then attempts to find a Dockerfile for these tools (can be explicitly
    specified with ``--dockerfile`` but by default it will check the tool's
    directory and the current directory as well).

    This command will then build and tag the image so it is ready to be tested
    and published. The docker_shell command be used to test out the built
    image.

    \b
        % planemo docker_build bowtie2.xml # asssumes Dockerfile in same dir
        % planemo docker_shell --from_tag bowtie2.xml

    This can optionally also cache the images.
    """
    dockerfiles.dockerfile_build(path=path, dockerfile=dockerfile, error=error, **kwds)
