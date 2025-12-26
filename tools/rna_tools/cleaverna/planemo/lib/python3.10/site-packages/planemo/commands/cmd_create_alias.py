"""Module describing the planemo ``create_alias`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.io import info

try:
    import namesgenerator
except ImportError:
    namesgenerator = None


@click.command("create_alias")
@click.argument(
    "obj",
    metavar="OBJ",
    type=click.STRING,
)
@options.alias_option()
@options.profile_option(required=True)
@command_function
def cli(ctx, alias, obj, profile, **kwds):
    """
    Add an alias for a path or a workflow or dataset ID. Aliases are associated with a particular planemo profile.
    """
    if not alias:
        if not namesgenerator:
            raise ImportError(
                "Random generation of aliases requires installation of the namesgenerator package."
                "Either install this, or specify the alias name with --alias."
            )
        alias = namesgenerator.get_random_name()

    exit_code = profiles.create_alias(ctx, alias, obj, profile)
    info(f"Alias {alias} created.")
    ctx.exit(exit_code)
    return
