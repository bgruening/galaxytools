"""Module describing the planemo ``create_gist`` command."""

import click

from planemo import github_util
from planemo.cli import command_function
from planemo.io import info

target_path = click.Path(
    file_okay=True,
    dir_okay=False,
    resolve_path=True,
)


@click.command("create_gist")
@click.argument(
    "path",
    metavar="FILE_PATH",
    type=target_path,
)
@click.option(
    "--link_type", type=click.Choice(["raw", "html"]), default="raw", help=("Link type to generate for the file.")
)
@command_function
def cli(ctx, path, **kwds):
    """Upload file to GitHub as a sharable gist."""
    file_url = github_util.publish_as_gist_file(ctx, path)
    if kwds.get("link_type") == "raw":
        share_url = file_url
    else:
        share_url = "http://htmlpreview.github.io/?%s" % file_url
    info("File published to Github Gist - share with %s" % share_url)
