"""Ensure best-practice biocontainer registered for this tool."""

from galaxy.tool_util.deps.container_resolvers.mulled import targets_to_mulled_name
from galaxy.tool_util.deps.mulled.util import build_target

from planemo.conda import tool_source_conda_targets

MESSAGE_WARN_NO_REQUIREMENTS = "No valid package requirement tags found to infer BioContainer from."
MESSAGE_WARN_NO_CONTAINER = "Failed to find a BioContainer registered for these requirements."
MESSAGE_INFO_FOUND_BIOCONTAINER = "BioContainer best-practice container found [%s]."

lint_tool_types = ["*"]


def lint_biocontainer_registered(tool_source, lint_ctx):
    conda_targets = tool_source_conda_targets(tool_source)
    if not conda_targets:
        lint_ctx.warn(MESSAGE_WARN_NO_REQUIREMENTS)
        return

    mulled_targets = [build_target(c.package, c.version) for c in conda_targets]
    name = mulled_container_name("biocontainers", mulled_targets)
    if name:
        lint_ctx.info(MESSAGE_INFO_FOUND_BIOCONTAINER % name)
    else:
        lint_ctx.warn(MESSAGE_WARN_NO_CONTAINER)


def mulled_container_name(namespace, targets):
    name = targets_to_mulled_name(targets=targets, hash_func="v2", namespace=namespace)
    if name:
        return f"quay.io/{namespace}/{name}"
