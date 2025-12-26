"""Module describing the planemo ``ci_find_repos`` command."""

import click

from planemo import options
from planemo.ci import (
    filter_paths,
    print_path_list,
)
from planemo.cli import command_function
from planemo.shed import find_raw_repositories


@click.command("ci_find_repos")
@options.shed_project_arg()
@options.ci_find_options()
@command_function
def cli(ctx, paths, **kwds):
    """Find all shed repositories in one or more directories.

    Currently, a repository is considered any directory with a .shed.yml
    or .dockstore.yml file.
    """
    kwds["recursive"] = True
    kwds["fail_fast"] = True
    repos = find_raw_repositories(ctx, paths, **kwds)
    # Since fail_fast is True, all repos are actual raw repo objects and
    # not exceptions.
    raw_paths = [r.path for r in repos]
    paths = filter_paths(ctx, raw_paths, path_type="repo", **kwds)
    print_path_list(paths, **kwds)
