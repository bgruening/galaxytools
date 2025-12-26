"""Module describing the planemo ``upload_data`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.galaxy.activity import stage_in
from planemo.galaxy.workflows import rewrite_job_file
from planemo.io import info
from planemo.runnable_resolve import for_runnable_identifier


@click.command("upload_data")
@options.required_runnable_arg()
@options.required_job_arg()
@options.required_new_job_arg()
@options.galaxy_run_options()
@options.galaxy_config_options()
@options.run_history_tags_option()
@command_function
def cli(ctx, runnable_identifier, job_path, new_job_path, **kwds):
    """Planemo command for uploading data to an external Galaxy server.

    \b
        % planemo upload_data wf.ga wf-job.yml new-wf-job.yml --profile profile

    Running this subcommand requires a workflow file or identifier
    and a job file, just as for ``planemo run``. In addition, a third
    argument is required, for the location of a new job file which
    Planemo will write. The data will be uploaded to the specified
    external Galaxy server and the job file will be recreated at the
    specified location, with all instances of ``path`` or ``location``
    for input datasets and collections replaced by ``galaxy_id``. The
    new job file can then be used to run the workflow separately from
    the already completed data upload.
    """
    kwds["engine"] = "external_galaxy"  # to only upload data to a local Galaxy does not make sense
    runnable = for_runnable_identifier(ctx, runnable_identifier, kwds)
    with engine_context(ctx, **kwds) as engine:
        with engine.ensure_runnables_served([runnable]) as config:
            job, _ = stage_in(ctx, runnable, config, job_path, **kwds)

    rewrite_job_file(job_path, new_job_path, job)
    info(f"Files uploaded and new job file written to {new_job_path}")

    return 0
