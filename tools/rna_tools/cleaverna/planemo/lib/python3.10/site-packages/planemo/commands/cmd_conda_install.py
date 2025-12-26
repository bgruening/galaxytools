"""Module describing the planemo ``conda_install`` command."""

import click
from galaxy.tool_util.deps import conda_util

from planemo import options
from planemo.cli import command_function
from planemo.conda import (
    build_conda_context,
    collect_conda_targets,
)
from planemo.io import coalesce_return_codes


@click.command("conda_install")
@options.optional_tools_or_packages_arg(multiple=True)
@options.recursive_option()
@options.conda_target_options()
@options.conda_global_option()
@options.conda_auto_init_option()
@command_function
def cli(ctx, paths, **kwds):
    """Install conda packages for tool requirements."""
    conda_context = build_conda_context(ctx, handle_auto_init=True, **kwds)
    return_codes = []
    for conda_target in collect_conda_targets(ctx, paths, recursive=kwds["recursive"]):
        ctx.log("Install conda target %s" % conda_target)
        return_code = conda_util.install_conda_target(
            conda_target, conda_context=conda_context, skip_environment=kwds.get("global", False)
        )
        return_codes.append(return_code)
    exit_code = coalesce_return_codes(return_codes, assert_at_least_one=True)
    ctx.exit(exit_code)
