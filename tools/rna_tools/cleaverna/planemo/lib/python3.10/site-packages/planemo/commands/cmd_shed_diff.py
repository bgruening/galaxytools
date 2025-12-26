"""Module describing the planemo ``shed_diff`` command."""

import shutil
import sys
import tempfile
from xml.sax.saxutils import escape

import click

from planemo import (
    options,
    shed,
)
from planemo.cli import command_function
from planemo.io import captured_io_for_xunit
from planemo.reports.xunit_handler import handle_report_xunit_kwd


@click.command("shed_diff")
@options.shed_read_options()
@click.option(
    "-o",
    "--output",
    type=click.Path(file_okay=True, resolve_path=True),
    help="Send diff output to specified file.",
    default=None,
)
@click.option(
    "--shed_target_source",
    help="Source Tool Shed to diff against (will ignore local project info"
    " specified). To compare the main Tool Shed against the test, set"
    " this to testtoolshed.",
    default=None,
)
@click.option(
    "--raw",
    is_flag=True,
    help="Do not attempt smart diff of XML to filter out attributes populated by the Tool Shed.",
)
@options.report_xunit()
@command_function
def cli(ctx, paths, **kwds):
    """diff between local repository and Tool Shed.

    By default, this will produce a diff between this repository and what
    would be uploaded to the Tool Shed with the `shed_upload` command - but
    this command can be made to compare other combinations of repositories.
    Here are some examples

    \b
        $ # diff for this repository and the main Tool Shed
        $ planemo shed_diff
        $ # diff for this repository and the test Tool Shed
        $ planemo shed_diff --shed_target testtoolshed
        $ # diff for the test Tool Shed and main Tool Shed
        $ planemo shed_diff --shed_target_source testtoolshed
        $ # diff for two an explicitly specified repositories (ignores
        $ # current project's shed YAML file.)
        $ planemo shed_diff --owner peterjc --name blast_rbh
            --shed_target_source testtoolshed

    This command will return an exit code of:

    - 0 if there are no detected differences.
    - 1 if there are differences.
    - 2 if the target repository doesn't exist.
    - >200 if there are errors attempting to perform a diff.

    **Warning:** ``shed_diff`` doesn't inspect repository metadata, this
    difference applies only to the file contents of files that would actually be
    uploaded to the repository.
    """

    # In a little bit of cheating, we're defining this variable here to collect
    # a "report" on our shed_diff
    collected_data = {
        "results": {
            "total": 0,
            "errors": 0,
            "failures": 0,
            "skips": 0,
        },
        "suitename": "shed_diff",
        "tests": [],
    }

    def diff(realized_repository):
        # We create a temporary redirection from kwds's
        # output to our tempfile. This lets us capture the
        # diff and redirect it to their requested location as
        # well as to the XUnit report.
        diff_output = tempfile.NamedTemporaryFile(mode="r")
        user_requested_output = kwds.get("output", None)
        # Replace their output handle with ours
        kwds["output"] = diff_output.name

        captured_io = {}
        with captured_io_for_xunit(kwds, captured_io):
            result = shed.diff_repo(ctx, realized_repository, **kwds)

        # May be extraneous but just want to ensure entire file is written
        # before a copy is made.
        diff_output.flush()
        # Redirect a copy to user_requested_output if they did:
        if user_requested_output is not None:
            shutil.copy(diff_output.name, user_requested_output)
        else:
            with open(diff_output.name) as f:
                sys.stdout.write(f.read())
        # Rewind to the start of the file and read it in its entirety
        diff_output.seek(0)
        diff_output_contents = diff_output.read()
        diff_output.close()

        # Collect data about what happened
        collected_data["results"]["total"] += 1
        xunit_case = {
            "name": "shed-diff",
            "classname": realized_repository.name,
            "time": captured_io["time"],
            "stdout": captured_io["stdout"],
            "stderr": captured_io["stderr"],
        }
        if result >= 200:
            collected_data["results"]["errors"] += 1
            xunit_case.update(
                {
                    "errorType": "DiffError",
                    "errorMessage": "Error diffing repositories",
                    "errorContent": escape(diff_output_contents),
                    "time": captured_io["time"],
                }
            )
        elif result > 2:
            collected_data["results"]["failures"] += 1
            xunit_case.update(
                {
                    "errorType": "PlanemoDiffError",
                    "errorMessage": "Planemo error diffing repositories",
                    "errorContent": escape(diff_output_contents),
                }
            )
        elif result == 2:
            collected_data["results"]["failures"] += 1
            xunit_case.update(
                {
                    "errorType": "RepoDoesNotExist",
                    "errorMessage": "Target Repository does not exist",
                    "errorContent": escape(diff_output_contents),
                }
            )
        elif result == 1:
            collected_data["results"]["failures"] += 1
            xunit_case.update(
                {
                    "errorType": "Different",
                    "errorMessage": "Repository is different",
                    "errorContent": escape(diff_output_contents),
                }
            )

        # Append our xunit test case
        collected_data["tests"].append(xunit_case)

        return result

    exit_code = shed.for_each_repository(ctx, diff, paths, **kwds)

    handle_report_xunit_kwd(kwds, collected_data)

    sys.exit(exit_code)
