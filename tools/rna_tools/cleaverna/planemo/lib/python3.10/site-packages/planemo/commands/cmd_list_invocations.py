"""Module describing the planemo ``list_invocations`` command."""

import json

import click

from planemo import options
from planemo.cli import command_function
from planemo.galaxy import profiles
from planemo.galaxy.api import get_invocations
from planemo.galaxy.workflows import remote_runnable_to_workflow_id
from planemo.io import (
    error,
    info,
)
from planemo.runnable_resolve import for_runnable_identifier

try:
    from tabulate import tabulate
except ImportError:
    tabulate = None  # type: ignore


@click.command("list_invocations")
@click.argument(
    "workflow_identifier",
    type=click.STRING,
)
@options.profile_option(required=True)
@command_function
def cli(ctx, workflow_identifier, **kwds):
    """
    Get a list of invocations for a particular workflow ID or alias.
    """
    info(f"Looking for invocations for workflow {workflow_identifier}...")
    runnable = for_runnable_identifier(ctx, workflow_identifier, kwds)
    profile = profiles.ensure_profile(ctx, kwds.get("profile"))
    assert runnable.is_remote_workflow_uri
    workflow_id = remote_runnable_to_workflow_id(runnable)

    invocations = get_invocations(
        url=profile["galaxy_url"],
        key=profile["galaxy_admin_key"] or profile["galaxy_user_key"],
        workflow_id=workflow_id,
    )
    if tabulate is not None:
        state_colors = {
            "ok": "\033[92m",  # green
            "running": "\033[93m",  # yellow
            "error": "\033[91m",  # red
            "paused": "\033[96m",  # cyan
            "deleted": "\033[95m",  # magenta
            "deleted_new": "\033[95m",  # magenta
            "new": "\033[96m",  # cyan
            "queued": "\033[93m",  # yellow
        }
        print(
            tabulate(
                {
                    "Invocation ID": invocations.keys(),
                    "Jobs status": [
                        ", ".join([f"{state_colors[k]}{v} jobs {k}\033[0m" for k, v in inv["states"].items()])
                        for inv in invocations.values()
                    ],
                    "Invocation report URL": [
                        "{}/workflows/invocations/report?id={}".format(profile["galaxy_url"].strip("/"), inv_id)
                        for inv_id in invocations
                    ],
                    "History URL": [
                        "{}/histories/view?id={}".format(
                            profile["galaxy_url"].strip("/"), invocations[inv_id]["history_id"]
                        )
                        for inv_id in invocations
                    ],
                },
                headers="keys",
            )
        )
    else:
        error("The tabulate package is not installed, invocations could not be listed correctly.")
        print(json.dumps(invocations, indent=4, sort_keys=True))
    info(f"{len(invocations)} invocations found.")
    return
