"""Ensure requirements are matched in best practice conda channels."""

from planemo.conda import (
    BEST_PRACTICE_CHANNELS,
    best_practice_search,
    tool_source_conda_targets,
)

lint_tool_types = ["*"]


def lint_requirements_in_conda(tool_source, lint_ctx):
    """Check requirements of tool source against best practice Conda channels."""
    conda_targets = tool_source_conda_targets(tool_source)
    if not conda_targets:
        lint_ctx.warn("No valid package requirement tags found to check against Conda.")
        return

    for conda_target in conda_targets:
        (best_hit, exact) = best_practice_search(conda_target)
        conda_target_str = conda_target.package
        if conda_target.version:
            conda_target_str += "@%s" % (conda_target.version)
        if best_hit and exact:
            template = "Requirement [%s] matches target in best practice Conda channel [%s]."
            message = template % (conda_target_str, best_hit.get("channel"))
            lint_ctx.info(message)
        elif best_hit:
            template = (
                "Requirement [%s] doesn't exactly match available version [%s] in best practice Conda channel [%s]."
            )
            message = template % (conda_target_str, best_hit["version"], best_hit.get("channel"))
            lint_ctx.warn(message)
        else:
            template = "Requirement [%s] doesn't match any recipe in a best practice conda channel [%s]."
            message = template % (conda_target_str, BEST_PRACTICE_CHANNELS)
            lint_ctx.warn(message)
