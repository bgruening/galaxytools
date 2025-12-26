"""Module describing the planemo ``workflow_upload`` command."""

from collections import defaultdict
from pathlib import Path

import click
from gxformat2.yaml import ordered_load_path

from planemo import options
from planemo.cli import command_function
from planemo.github_util import create_release
from planemo.workflow_lint import find_workflow_descriptions


@click.command("workflow_upload")
@options.github_namespace()
@options.github_branch()
@options.dry_run()
@options.optional_tools_or_packages_arg(multiple=True)
@command_function
def cli(ctx, paths, namespace, dry_run, github_branch, **kwds):
    """Upload workflows to github organization."""
    owner = namespace
    for path in paths:
        path = Path(path).absolute()
        if path.is_dir():
            repo = path.name
        else:
            repo = path.parent.name
            path = path.parent
        versions = defaultdict(list)
        for workflow_file in find_workflow_descriptions(path):
            workflow = ordered_load_path(workflow_file)
            version = workflow.get("release")
            if not version:
                raise Exception(f"Must set a release version in workflow file '{workflow_file}'")
            versions[version].append(workflow_file)
            if len(versions) > 1:
                msg = ""
                for version, paths in versions.items():
                    msg = "{}version: {}\npaths: {}".format(msg, version, "\n".join(paths))
                raise Exception(f"All workflows in repository must have same version.\n{msg}")
        if versions:
            create_release(
                ctx,
                from_dir=path,
                target_dir=repo,
                owner=owner,
                repo=repo,
                version=version,
                dry_run=dry_run,
                branch=github_branch,
            )
