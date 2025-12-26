"""Module describing the planemo ``shed_lint`` command."""

import click

from planemo import (
    options,
    shed,
    shed_lint,
)
from planemo.cli import (
    command_function,
    PlanemoCliContext,
)


@click.command("shed_lint")
@options.shed_realization_options()
@options.report_level_option()
@options.fail_level_option()
@options.skip_options()
@options.click.option(
    "--tools", is_flag=True, default=False, help=("Lint tools discovered in the process of linting repositories.")
)
@options.click.option(
    "--ensure_metadata",
    is_flag=True,
    default=False,
    help=(
        "Ensure .shed.yml files contain enough metadata for each repository to allow automated creation and/or updates."
    ),
)
@click.option(
    "--urls",
    is_flag=True,
    default=False,
    help="Check validity of URLs in XML files",
)
@options.lint_biocontainers_option()
# @click.option(
#     "--verify",
#     is_flag=True,
#     help="If an sha256sum is available, download the entire file AND validate it.",
#     default=False,
# )
@command_function
def cli(ctx: PlanemoCliContext, paths, **kwds):
    """Check Tool Shed repository for common issues.

    With the ``--tools`` flag, this command lints actual Galaxy tools
    in addition to tool shed artifacts.

    With the ``--urls`` flag, this command searches for
    ``<package>$URL</package>`` and download actions which specify URLs. Each
    of those are accessed individually. By default, this tool requests the
    first hundred or so bytes of each listed URL and validates that a 200 OK
    was received. In tool XML files, the ``--urls`` option checks through the
    help text for mentioned URLs and checks those.
    """

    def lint(realized_repository):
        return shed_lint.lint_repository(ctx, realized_repository, **kwds)

    kwds["fail_on_missing"] = False
    exit_code = shed.for_each_repository(ctx, lint, paths, **kwds)
    ctx.exit(exit_code)
