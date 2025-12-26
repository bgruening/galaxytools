"""Module describing the planemo ``list_repos`` command."""

import click

from planemo import options
from planemo.ci import print_as_yaml
from planemo.cli import command_function
from planemo.shed import (
    _realize_effective_repositories,
    find_raw_repositories,
)


@click.command("list_repos")
@options.shed_project_arg()
@options.ci_find_options()
@command_function
def cli(ctx, paths, **kwds):
    """Find all shed repositories in one or more directories and output as yaml.

    Currently, a shed repository is considered a directory with a .shed.yml
    file.
    """
    kwds["recursive"] = True
    kwds["fail_fast"] = True
    repos = find_raw_repositories(ctx, paths, **kwds)
    # Since fail_fast is True, all repos are actual raw repo objects and
    # not exceptions.
    raw_paths = [r.path for r in repos]
    realized_repos = []
    for repo_gen in (_realize_effective_repositories(r, path=p) for r, p in zip(repos, raw_paths)):
        for repo in repo_gen:
            realized_repos.append({"name": repo.config["name"], "owner": repo.config["owner"]})
    print_as_yaml(realized_repos, **kwds)
