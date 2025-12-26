"""Module describing the planemo ``workflow_test_check`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine.factory import engine_context
from planemo.galaxy.activity import invocation_to_run_response
from planemo.galaxy.test.actions import handle_reports_and_summary
from planemo.galaxy.workflows import GALAXY_WORKFLOW_INSTANCE_PREFIX
from planemo.runnable import definition_to_test_case
from planemo.runnable_resolve import for_runnable_identifier
from planemo.test.results import StructuredData


@click.command("workflow_test_on_invocation")
@options.optional_tools_arg(multiple=False, allow_uris=False, metavar="TEST.YML")
@options.invocation_target_options()
@options.test_index_option()
@options.test_output_options()
@command_function
def cli(ctx, path, invocation_id, test_index, **kwds):
    """Run defined tests against existing workflow invocation."""
    with engine_context(ctx, engine="external_galaxy", **kwds) as engine, engine.ensure_runnables_served([]) as config:
        user_gi = config.user_gi
        invocation = user_gi.invocations.show_invocation(invocation_id)
        runnable = for_runnable_identifier(
            ctx, f"{GALAXY_WORKFLOW_INSTANCE_PREFIX}{invocation['workflow_id']}?runnable_path={path}", kwds
        )
        test_cases = definition_to_test_case(path, runnable)
        assert (
            len(test_cases) >= test_index
        ), f"Selected test case {test_index}, but only found {len(test_cases)} test case(s)."
        test_case = test_cases[test_index - 1]
        # Hardcode fail_fast, no need to expose the option to the user IMO.
        run_response = invocation_to_run_response(
            ctx, user_gi=config.user_gi, runnable=runnable, invocation=invocation, fail_fast=True
        )
        structured_data = test_case.structured_test_data(run_response)
        test_data = {
            "version": "0.1",
            "tests": [structured_data],
        }
        structured_results = StructuredData(data=test_data)
        structured_results.calculate_summary_data()
    return_value = handle_reports_and_summary(ctx, structured_results.structured_data, kwds=kwds)
    ctx.exit(return_value)
