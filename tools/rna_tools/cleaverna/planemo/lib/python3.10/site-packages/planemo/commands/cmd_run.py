"""Module describing the planemo ``run`` command."""

import json

import click
from galaxy.util import unicodify

from planemo import options
from planemo.cli import command_function
from planemo.engine import engine_context
from planemo.galaxy.test import handle_reports_and_summary
from planemo.io import (
    info,
    warn,
)
from planemo.runnable import RunnableType
from planemo.runnable_resolve import for_runnable_identifier
from planemo.test.results import StructuredData


@click.command("run")
@options.required_runnable_arg()
@options.required_job_arg()
@options.galaxy_run_options()
@options.galaxy_config_options()
@options.enable_cwl_option()
@options.galaxy_cwl_root_option()
@options.run_history_tags_option()
@options.run_output_directory_option()
@options.run_output_json_option()
@options.run_download_outputs_option()
@options.run_export_option()
@options.invocation_export_format_arg()
@options.engine_options()
@options.test_options()
@command_function
def cli(ctx, runnable_identifier, job_path, **kwds):
    """Planemo command for running tools and jobs.

    \b
        % planemo run cat1-tool.cwl cat-job.json
    """
    runnable = for_runnable_identifier(ctx, runnable_identifier, kwds)
    is_cwl = runnable.type.is_cwl_artifact
    if kwds.get("export_invocation") and not runnable.type == RunnableType.galaxy_workflow:
        raise click.UsageError(
            "Exporting invocation is only supported for Galaxy workflows, "
            "but the provided runnable is of type: %s" % runnable.type
        )
    kwds["cwl"] = is_cwl
    kwds["execution_type"] = "Run"
    if kwds.get("engine", None) is None:
        if is_cwl:
            kwds["engine"] = "cwltool"
        elif kwds.get("galaxy_url", None):
            kwds["engine"] = "external_galaxy"
        else:
            kwds["engine"] = "galaxy"
    with engine_context(ctx, **kwds) as engine:
        run_result = engine.run([runnable], [job_path])[0]

    if not run_result.was_successful:
        warn("Run failed [%s]" % unicodify(run_result))
    elif kwds.get("no_wait"):
        info("Run successfully executed - exiting without waiting for results.")
    else:
        output_json = kwds.get("output_json", None)
        outputs_dict = run_result.outputs_dict
        if output_json:
            with open(output_json, "w") as f:
                json.dump(outputs_dict, f, ensure_ascii=False)
        info("Run completed successfully.")

    report_data = StructuredData(data={"tests": [run_result.structured_data()], "version": "0.1"})
    report_data.calculate_summary_data()
    return_value = handle_reports_and_summary(ctx, report_data.structured_data, kwds=kwds)
    ctx.exit(return_value)
