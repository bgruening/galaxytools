"""Module describing the planemo ``test`` command."""

import click

from planemo import options
from planemo.cli import command_function
from planemo.engine.test import test_runnables
from planemo.runnable import RunnableType
from planemo.runnable_resolve import for_runnable_identifiers


@click.command("test")
@options.optional_tools_arg(multiple=True, allow_uris=True)
@click.option(
    "--failed",
    is_flag=True,
    help="Re-run only failed tests. This command will read "
    "tool_test_output.json to determine which tests failed so this "
    "file must have been produced with the same set of tool ids "
    "previously.",
    default=False,
)
@click.option(
    "--polling_backoff",
    type=int,
    help="Poll resources with an increasing interval between requests. "
    "Useful when testing against remote and/or production "
    "instances to limit generated traffic.",
    default="0",
)
@options.galaxy_target_options()
@options.galaxy_config_options()
@options.test_options()
@options.engine_options()
@command_function
def cli(ctx, uris, **kwds):
    """Run specified tool or workflow tests within Galaxy.

    All referenced tools (by default all the tools in the current working
    directory) will be tested and the results quickly summarized.

    To run these tests planemo needs a Galaxy instance to utilize, planemo
    will search parent directories to see if any is a Galaxy instance
    - but one can pick the Galaxy instance to use with the --galaxy_root
    option or force planemo to download a disposable instance with the
    ``--install_galaxy`` flag.

    In addition to to quick summary printed to the console - various detailed
    output summaries can be configured. ``tool_test_output.html`` (settable
    via ``--test_output``) will contain a human consumable HTML report
    describing the test run. A JSON file (settable via ``--test_output_json``
    and defaulting to ``tool_test_output.json``) will also be created. These
    files can can be disabled by passing in empty arguments or globally by
    setting the values ``default_test_output`` and/or
    ``default_test_output_json`` in ``~/.planemo.yml`` to ``null``. For
    continuous integration testing a xUnit-style report can be configured using
    the ``--test_output_xunit``.

    planemo uses temporarily generated config files and environment variables
    to attempt to shield this execution of Galaxy from manually launched runs
    against that same Galaxy root - but this may not be bullet proof yet so
    please careful and do not try this against production Galaxy instances.
    """
    runnables = for_runnable_identifiers(ctx, uris, kwds)

    # pick a default engine type if needed
    is_cwl = all(r.type in {RunnableType.cwl_tool, RunnableType.cwl_workflow} for r in runnables)
    if kwds.get("engine", None) is None:
        if is_cwl:
            kwds["engine"] = "cwltool"
        elif kwds.get("galaxy_url", None):
            kwds["engine"] = "external_galaxy"
        else:
            kwds["engine"] = "galaxy"

    return_value = test_runnables(ctx, runnables, original_paths=uris, **kwds)

    ctx.exit(return_value)
