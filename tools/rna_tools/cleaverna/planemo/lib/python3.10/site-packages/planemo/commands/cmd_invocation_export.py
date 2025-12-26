"""Module describing the planemo ``invocation_export`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.api import export_invocation_as_archive
from planemo.io import info


@click.command("invocation_export")
@options.required_invocation_id_arg()
@options.invocation_export_format_arg()
@options.ci_output_option()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@options.profile_option()
@command_function
def cli(ctx, invocation_id, export_format, output, galaxy_url, galaxy_user_key, **kwds):
    """Planemo command for exporting an invocation as an archive.
    \b
        % planemo invocation_export ID --profile my_profile --output invocation.rocrate.zip
    """

    export_invocation_as_archive(
        invocation_id=invocation_id,
        export_format=export_format,
        output=output,
        url=galaxy_url,
        key=galaxy_user_key,
    )
    info(f"Exported invocation {invocation_id} to {output}, format: {export_format}")
