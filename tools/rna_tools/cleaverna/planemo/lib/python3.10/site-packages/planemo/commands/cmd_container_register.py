"""Module describing the planemo ``container_register`` command."""

import os
import string
from typing import List

import click
from galaxy.tool_util.deps.container_resolvers.mulled import targets_to_mulled_name
from galaxy.tool_util.deps.mulled.mulled_build import (
    base_image_for_targets,
    DEFAULT_BASE_IMAGE,
)
from galaxy.tool_util.deps.mulled.util import (
    conda_build_target_str,
    CondaTarget,
    v2_image_name,
)

from planemo import options
from planemo.cli import (
    command_function,
    PlanemoCliContext,
)
from planemo.conda import (
    best_practice_search,
    build_conda_context,
    collect_conda_target_lists_and_tool_paths,
)
from planemo.git import (
    add,
    branch,
    commit,
    push,
)
from planemo.github_util import (
    clone_fork_branch,
    DEFAULT_REMOTE_NAME,
    get_repository_object,
    pull_request,
)
from planemo.mulled import conda_to_mulled_targets

REGISTRY_TARGET_NAME = "multi-package-containers"
REGISTRY_TARGET_PATH = "combinations"
REGISTRY_REPOSITORY = "BioContainers/multi-package-containers"
DEFAULT_MESSAGE = "Add container $hash.\n**Hash**: $hash\n\n**Packages**:\n$packages\nBase Image:$base_image\n\n"
DEFAULT_MESSAGE += "**For** :\n$tools\n\nGenerated with Planemo."
CONTENTS = "#targets\tbase_image\timage_build\n$targets\t$base_image\t$image_build\n"
BIOCONTAINERS_PLATFORM = "linux-64"


@click.command("container_register")
@options.optional_tools_arg(multiple=True)
@options.recursive_option()
@options.mulled_namespace_option()
@options.conda_target_options()
@click.option(
    "output_directory",
    "--output_directory",
    type=click.Path(
        file_okay=False,
        dir_okay=True,
        resolve_path=True,
    ),
    default=None,
    help=("Container registration directory (defaults to ~/.planemo/multi-package-containers."),
)
@click.option(
    "-m",
    "--message",
    default=DEFAULT_MESSAGE,
    help="Commit and pull request message template for registration interactions.",
)
@click.option(
    "--pull_request/--no_pull_request",
    is_flag=True,
    default=True,
    help="Fork and create a pull request against %s for these changes." % REGISTRY_REPOSITORY,
)
@click.option(
    "--force_push/--no_force_push",
    is_flag=True,
    default=False,
    help="Force push branch for pull request in case it already exists.",
)
@command_function
def cli(ctx: "PlanemoCliContext", paths, **kwds) -> None:
    """Register multi-requirement containers as needed.

    BioContainers publishes all Bioconda packages automatically as individual
    container images. These however are not enough for tools with multiple
    best-practice requirements. Such requirements should be recorded and published
    so that a container can be created and registered for these tools.
    """
    registry_target = RegistryTarget(ctx, **kwds)
    conda_context = build_conda_context(ctx, **kwds)

    combinations_added = 0
    conda_targets_list, tool_paths_list = collect_conda_target_lists_and_tool_paths(
        ctx, paths, recursive=kwds["recursive"]
    )
    for conda_targets, tool_paths in zip(conda_targets_list, tool_paths_list):
        ctx.vlog("Handling conda_targets [%s]" % conda_targets)
        mulled_targets = conda_to_mulled_targets(conda_targets)
        mulled_targets_str = "- " + "\n- ".join(map(conda_build_target_str, mulled_targets))

        if len(mulled_targets) < 1:
            ctx.log("Skipping registration, no targets discovered.")
            continue

        name = v2_image_name(mulled_targets)
        tag = "0"
        name_and_tag = f"{name}-{tag}"
        target_filename = os.path.join(registry_target.output_directory, "%s.tsv" % name_and_tag)
        if os.path.exists(target_filename):
            ctx.log("Target file '%s' already exists, skipping" % target_filename)
            continue

        if targets_to_mulled_name(mulled_targets, hash_func="v2", namespace=kwds["mulled_namespace"]):
            ctx.vlog("quay repository already exists, skipping")
            continue

        if registry_target.has_pull_request_for(name):
            ctx.vlog("Found matching open pull request for [%s], skipping" % name)
            continue

        best_practice_requirements = True
        for conda_target in conda_targets:
            best_hit, exact = best_practice_search(
                conda_target, conda_context=conda_context, platform=BIOCONTAINERS_PLATFORM
            )
            if not best_hit:
                ctx.log("Target [%s] is not available in best practice channels - skipping" % conda_target)
                best_practice_requirements = False
            if not exact:
                ctx.log(
                    "Target version [%s] for package [%s] is not available in best practice channels - please specify the full version",
                    conda_target.version,
                    conda_target.package,
                )

        if not best_practice_requirements:
            continue

        base_image = DEFAULT_BASE_IMAGE
        for conda_target in conda_targets:
            base_image = base_image_for_targets([conda_target], conda_context=conda_context)
            if base_image != DEFAULT_BASE_IMAGE:
                ctx.log(f"{conda_target} requires '{base_image}' as base image")
                break

        registry_target.write_targets(ctx, target_filename, mulled_targets, tag, base_image)
        tools_str = "\n".join(map(lambda p: "- " + os.path.basename(p), tool_paths))
        registry_target.handle_pull_request(
            ctx, name, target_filename, mulled_targets_str, tools_str, base_image, **kwds
        )
        combinations_added += 1


class RegistryTarget:
    """Abstraction around mulled container registry (both directory and Github repo)."""

    def __init__(self, ctx: "PlanemoCliContext", **kwds):
        output_directory = kwds["output_directory"]
        pr_titles = []
        target_repository = None
        self.remote_name = DEFAULT_REMOTE_NAME
        do_pull_request = kwds.get("pull_request", True)
        if output_directory is None:
            target_repository = os.path.join(ctx.workspace, REGISTRY_TARGET_NAME)
            output_directory = os.path.join(target_repository, REGISTRY_TARGET_PATH)
            self.remote_name = (
                clone_fork_branch(
                    ctx,
                    "https://github.com/%s" % REGISTRY_REPOSITORY,
                    target_repository,
                    fork=do_pull_request,
                )
                or self.remote_name
            )
            pr_titles = [pr.title for pr in open_prs(ctx)]

        self.do_pull_request = do_pull_request
        self.pr_titles = pr_titles
        self.output_directory = output_directory
        self.target_repository = target_repository

    def has_pull_request_for(self, name: str) -> bool:
        has_pr = False
        if self.do_pull_request:
            if any([name in t for t in self.pr_titles]):
                has_pr = True

        return has_pr

    def handle_pull_request(
        self,
        ctx: "PlanemoCliContext",
        name: str,
        target_filename: str,
        packages_str: str,
        tools_str: str,
        base_image: str,
        **kwds,
    ) -> None:
        if self.do_pull_request:
            message = kwds["message"]
            message = string.Template(message).safe_substitute(
                {
                    "hash": name,
                    "packages": packages_str,
                    "tools": tools_str,
                    "base_image": base_image,
                }
            )
            branch_name = name.replace(":", "-")
            assert self.target_repository
            branch(ctx, self.target_repository, branch_name, from_branch="master")
            add(ctx, self.target_repository, target_filename)
            commit(ctx, self.target_repository, message=message)
            force_push = kwds.get("force_push", False)
            push(ctx, repo_path=self.target_repository, to=self.remote_name, branch=branch_name, force=force_push)
            pull_request(ctx, self.target_repository, message=message, repo=REGISTRY_REPOSITORY)

    def write_targets(
        self,
        ctx: "PlanemoCliContext",
        target_filename: str,
        mulled_targets: List[CondaTarget],
        tag: str,
        base_image: str,
    ) -> None:
        with open(target_filename, "w") as f:
            targets = to_target_str(mulled_targets)
            f.write(string.Template(CONTENTS).safe_substitute(targets=targets, base_image=base_image, image_build=tag))
            ctx.log(f"Wrote requirements [{targets}] to file [{target_filename}]")


def to_target_str(targets: List[CondaTarget]) -> str:
    target_strings = []
    for target in targets:
        if target.version:
            target_str = f"{target.package}={target.version}"
        else:
            target_str = target.package
        target_strings.append(target_str)
    return ",".join(target_strings)


def open_prs(ctx: "PlanemoCliContext") -> List:
    repo = get_repository_object(ctx, REGISTRY_REPOSITORY)
    prs = [pr for pr in repo.get_pulls()]
    return prs
