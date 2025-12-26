"""Module describing the planemo ``autoupdate`` command."""

import json
from typing import Callable

import click
import yaml
from galaxy.tool_util.parser.xml import XmlToolSource
from gxformat2.export import from_galaxy_native
from gxformat2.yaml import ordered_dump

from planemo import (
    autoupdate,
    options,
)
from planemo.cli import command_function
from planemo.config import planemo_option
from planemo.engine import (
    engine_context,
    is_galaxy_engine,
)
from planemo.engine.test import test_runnables
from planemo.exit_codes import (
    EXIT_CODE_GENERIC_FAILURE,
    EXIT_CODE_OK,
)
from planemo.io import (
    coalesce_return_codes,
    error,
    info,
)
from planemo.runnable import (
    for_paths,
    RunnableType,
)
from planemo.tools import (
    is_tool_load_error,
    yield_tool_sources_on_paths,
)


def dry_run_option() -> Callable:
    """Perform a dry run autoupdate without modifying the XML files"""
    return planemo_option(
        "--dry-run", is_flag=True, help="Perform a dry run autoupdate without modifying the XML files."
    )


def test_option() -> Callable:
    """Test updated XML files"""
    return planemo_option("--test", is_flag=True, help="Test updated XML files.")


def skiplist_option() -> Callable:
    """List of XML files to skip"""
    return planemo_option(
        "--skiplist",
        default=None,
        help="Skiplist file, containing a list of tools or workflows for which autoupdate should be skipped.",
    )


def skip_requirements_option() -> Callable:
    """List of requirements to skip"""
    return planemo_option(
        "--skip_requirements",
        default="python,r-base,perl",
        help="Comma-separated list of requirements which should be not be updated. Default is python,r-base,perl.",
    )


@click.command("autoupdate")  # noqa C901
@options.optional_tools_arg(multiple=True)
@dry_run_option()
@options.recursive_option()
@test_option()
@skiplist_option()
@skip_requirements_option()
@options.serve_engine_option()
@options.test_options()
@options.galaxy_target_options()
@options.galaxy_config_options()
@options.report_level_option()
@options.report_xunit()
@options.fail_level_option()
@options.galaxy_url_option()
@options.galaxy_user_key_option()
@options.galaxy_admin_key_option()
@command_function
def cli(ctx, paths, **kwds):  # noqa C901
    """Auto-update tool requirements by checking against Conda and updating if newer versions are available."""
    assert_tools = kwds.get("assert_tools", True)
    recursive = kwds.get("recursive", False)
    exit_codes = []
    modified_files = set()
    tools_to_skip = [line.rstrip() for line in open(kwds["skiplist"])] if kwds["skiplist"] else []
    runnables = for_paths(paths)

    if any(r.type in {RunnableType.galaxy_tool, RunnableType.directory} for r in runnables):
        # update Galaxy tools
        for tool_path, tool_xml in yield_tool_sources_on_paths(ctx, paths, recursive):
            if tool_path.split("/")[-1] in tools_to_skip:
                info("Skipping tool %s" % tool_path)
                continue
            info("Auto-updating tool %s" % tool_path)
            try:
                updated = autoupdate.autoupdate_tool(ctx, tool_path, modified_files=modified_files, **kwds)
                if updated:
                    modified_files.update(updated)
            except Exception as e:
                error(f"{tool_path} could not be updated - the following error was raised: {e.__str__()}")
            if handle_tool_load_error(tool_path, tool_xml):
                exit_codes.append(EXIT_CODE_GENERIC_FAILURE)
                continue
            else:
                exit_codes.append(EXIT_CODE_OK)

    workflows = [
        r for r in runnables if r.type == RunnableType.galaxy_workflow and r.path.split("/")[-1] not in tools_to_skip
    ]
    modified_workflows = []
    for workflow in workflows:
        tools_to_update = autoupdate.get_tools_to_update(ctx, workflow, tools_to_skip)
        if tools_to_update:
            modified_workflows.append(workflow)
            info("The following tools are outdated:")
            for tool, versions in tools_to_update.items():
                info(f"{tool}: {versions['current']} -> {versions['updated']}")
        else:
            info("No newer tool versions were found, so the workflow was not updated.")

    if modified_workflows and not kwds.get("dry_run"):
        assert is_galaxy_engine(**kwds)
        if kwds.get("engine") != "external_galaxy" or kwds.get("galaxy_admin_key"):
            kwds["install_most_recent_revision"] = True
            kwds["install_resolver_dependencies"] = False
            kwds["install_repository_dependencies"] = False
            kwds["shed_install"] = True

        with engine_context(ctx, **kwds) as galaxy_engine:
            with galaxy_engine.ensure_runnables_served(modified_workflows) as config:
                for workflow in modified_workflows:
                    info("Auto-updating workflow %s" % workflow.path)
                    updated_workflow = autoupdate.autoupdate_wf(ctx, config, workflow)

                    with open(workflow.path) as f:
                        original_workflow = yaml.load(f, Loader=yaml.SafeLoader)
                    edited_workflow = autoupdate.fix_workflow(original_workflow, updated_workflow)
                    with open(workflow.path, "w") as f:
                        if workflow.path.endswith(".ga"):
                            json.dump(edited_workflow, f, ensure_ascii=False, indent=4)
                        else:
                            gxformat2_wf = from_galaxy_native(edited_workflow)
                            if "release" in edited_workflow:
                                gxformat2_wf["release"] = edited_workflow["release"]
                            ordered_dump(gxformat2_wf, f)
                    if original_workflow.get("release"):
                        info(
                            f"The workflow release number has been updated from "
                            f"{original_workflow.get('release')} to {edited_workflow.get('release')}."
                        )

    if kwds["test"]:
        if not modified_files and not modified_workflows:
            info("No tools or workflows were updated, so no tests were run.")
        else:
            modified_paths = [
                path for path, tool_xml in yield_tool_sources_on_paths(ctx, paths, recursive) if path in modified_files
            ]
            info(f"Running tests for the following auto-updated tools: {', '.join(modified_paths)}")
            runnables = for_paths(modified_paths + modified_workflows)
            kwds["engine"] = "galaxy"
            return_value = test_runnables(ctx, runnables, original_paths=paths, **kwds)
            exit_codes.append(return_value)
    return coalesce_return_codes(exit_codes, assert_at_least_one=assert_tools)


def handle_tool_load_error(tool_path: str, tool_xml: XmlToolSource) -> bool:
    """Return True if tool_xml is tool load error (invalid XML), and
    print a helpful error message.
    """
    is_error = False
    if is_tool_load_error(tool_xml):
        info("Could not update %s due to malformed xml." % tool_path)
        is_error = True
    return is_error
