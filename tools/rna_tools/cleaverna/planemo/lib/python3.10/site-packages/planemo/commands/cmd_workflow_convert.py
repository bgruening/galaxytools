"""Module describing the planemo ``workflow_convert`` command."""

import json
import os

import click
from gxformat2 import from_galaxy_native

from planemo import options
from planemo.cli import command_function
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.io import write_file
from planemo.runnable import for_path


@click.command("workflow_convert")
@options.required_workflow_arg()
@options.force_option()
@options.workflow_output_artifact()
@options.galaxy_serve_options()
@command_function
def cli(ctx, workflow_identifier, output=None, force=False, **kwds):
    """Convert Format 2 workflows to native Galaxy workflows, and vice-versa."""
    assert is_galaxy_engine(**kwds)

    kwds["no_dependency_resolution"] = True

    if workflow_identifier.endswith(".ga"):
        if output is None:
            output = os.path.splitext(workflow_identifier)[0] + ".gxwf.yml"

        with open(workflow_identifier) as f:
            workflow_dict = json.load(f)
        format2_wrapper = from_galaxy_native(workflow_dict, json_wrapper=True)
        with open(output, "w") as f:
            f.write(format2_wrapper["yaml_content"])
    else:
        if output is None:
            output = os.path.splitext(workflow_identifier)[0] + ".ga"

        runnable = for_path(workflow_identifier)
        with engine_context(ctx, **kwds) as galaxy_engine:
            with galaxy_engine.ensure_runnables_served([runnable]) as config:
                workflow_id = config.workflow_id(workflow_identifier)
                output_dict = config.gi.workflows.export_workflow_dict(workflow_id)

                output_contents = json.dumps(output_dict, ensure_ascii=False, indent=4, sort_keys=True)
                write_file(output, output_contents, force=force)
