"""Module describing the planemo ``workflow_test_init`` command."""

import os

import click
import yaml

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.workflows import (
    get_workflow_from_invocation_id,
    job_template,
    new_workflow_associated_path,
    output_stubs_for_workflow,
)
from planemo.io import can_write_to_path


@click.command("workflow_test_init")
@options.required_workflow_arg()
@options.force_option()
@options.workflow_output_artifact()
@options.split_job_and_test()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@options.from_invocation()
@options.profile_option()
@command_function
def cli(ctx, workflow_identifier, output=None, split_test=False, **kwds):
    """Initialize a Galaxy workflow test description for supplied workflow.

    Be sure to lint your workflow with ``workflow_lint`` before calling this
    to ensure inputs and outputs comply with best practices that make workflow
    testing easier.
    """
    if kwds["from_invocation"]:
        if not os.path.isdir("test-data"):
            ctx.log("Creating test-data directory.")
            os.makedirs("test-data")
        path_basename = get_workflow_from_invocation_id(
            workflow_identifier, kwds["galaxy_url"], kwds["galaxy_user_key"]
        )
    else:
        path_basename = os.path.basename(workflow_identifier)
    job = job_template(workflow_identifier, **kwds)
    if output is None:
        output = new_workflow_associated_path(path_basename if kwds["from_invocation"] else workflow_identifier)
    job_output = new_workflow_associated_path(
        path_basename if kwds["from_invocation"] else workflow_identifier, suffix="job1"
    )
    if not can_write_to_path(output, **kwds):
        ctx.exit(1)

    test_description = [
        {
            "doc": "Test outline for %s" % path_basename,
            "job": job,
            "outputs": output_stubs_for_workflow(workflow_identifier, **kwds),
        }
    ]
    if split_test:
        job_output = new_workflow_associated_path(
            path_basename if kwds["from_invocation"] else workflow_identifier, suffix="job1"
        )
        if not can_write_to_path(job_output, **kwds):
            ctx.exit(1)

        test_description[0]["job"] = os.path.basename(job_output)
        with open(job_output, "w") as f_job:
            yaml.dump(job, f_job, sort_keys=False)
    with open(output, "w") as f:
        yaml.dump(test_description, f, sort_keys=False)
