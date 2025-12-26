"""Module describing the planemo ``list_alias`` command."""

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.io import info

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None  # type: ignore


@click.command("list_alias")
@options.profile_option(required=True)
@command_function
def cli(ctx, profile, **kwds):
    """
    List aliases for a path or a workflow or dataset ID. Aliases are associated with a particular planemo profile.
    """
    info("Looking for profiles...")
    aliases = profiles.list_alias(ctx, profile)
    if tabulate is not None:
        print(tabulate({"Alias": aliases.keys(), "Object": aliases.values()}, headers="keys"))
    else:
        print(json.dumps(aliases, indent=4, sort_keys=True))

    info(f"{len(aliases)} aliases were found for profile {profile}.")

    ctx.exit(0)
    return
