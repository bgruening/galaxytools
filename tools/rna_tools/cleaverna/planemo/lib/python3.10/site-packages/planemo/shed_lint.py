"""Logic related to linting shed repositories."""

import os
import xml.etree.ElementTree as ET
from typing import TYPE_CHECKING

import yaml
from bioblend import ConnectionError
from galaxy.tool_util.lint import lint_tool_source_with
from galaxy.tool_util.linters.help import rst_invalid
from galaxy.tool_util.parser.interface import ToolSource
from galaxy.tool_util.version import parse_version
from galaxy.util import unicodify

from planemo.io import info
from planemo.lint import (
    handle_lint_complete,
    lint_urls,
    lint_xsd,
    setup_lint,
)
from planemo.shed import (
    CURRENT_CATEGORIES,
    REPO_TYPE_SUITE,
    REPO_TYPE_TOOL_DEP,
    REPO_TYPE_UNRESTRICTED,
    validate_repo_name,
    validate_repo_owner,
)
from planemo.shed2tap import base
from planemo.shed.interface import tool_shed_instance
from planemo.tool_lint import (
    build_tool_lint_args,
    handle_tool_load_error,
)
from planemo.tools import yield_tool_sources
from planemo.xml import XSDS_PATH

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext
    from planemo.shed import RealizedRepository

TOOL_DEPENDENCIES_XSD = os.path.join(XSDS_PATH, "tool_dependencies.xsd")
REPO_DEPENDENCIES_XSD = os.path.join(XSDS_PATH, "repository_dependencies.xsd")


VALID_REPOSITORY_TYPES = [
    REPO_TYPE_UNRESTRICTED,
    REPO_TYPE_TOOL_DEP,
    REPO_TYPE_SUITE,
]

SHED_METADATA = [
    "description",
    "long_description",
    "remote_repository_url",
    "homepage_url",
    "categories",
]


def lint_repository(ctx: "PlanemoCliContext", realized_repository: "RealizedRepository", **kwds):
    """Lint a realized shed repository.

    See :mod:`planemo.shed` for details on constructing a realized
    repository data structure.
    """
    failed = False
    path = realized_repository.real_path
    info("Linting repository %s" % path)
    lint_args = build_tool_lint_args(ctx, **kwds)
    lint_args, lint_ctx = setup_lint(ctx, lint_args=lint_args, **kwds)
    lint_ctx.lint(
        "lint_expansion",
        lint_expansion,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_expected_files",
        lint_expected_files,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_tool_dependencies_xsd",
        lint_tool_dependencies_xsd,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_tool_dependencies_sha256sum",
        lint_tool_dependencies_sha256sum,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_tool_dependencies_actions",
        lint_tool_dependencies_actions,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_repository_dependencies",
        lint_repository_dependencies,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_shed_yaml",
        lint_shed_yaml,
        realized_repository,
    )
    lint_ctx.lint(
        "lint_readme",
        lint_readme,
        realized_repository,
    )
    if kwds["urls"]:
        lint_ctx.lint(
            "lint_urls",
            lint_tool_dependencies_urls,
            realized_repository,
        )
    if kwds["tools"]:
        tools_failed = lint_repository_tools(ctx, realized_repository, lint_ctx, lint_args)
        failed = failed or tools_failed

    lint_ctx.lint("lint_version_bumped", lint_shed_version, realized_repository)

    if kwds["ensure_metadata"]:
        lint_ctx.lint(
            "lint_shed_metadata",
            lint_shed_metadata,
            realized_repository,
        )
    return handle_lint_complete(lint_ctx, lint_args, failed=failed)


def lint_repository_tools(ctx: "PlanemoCliContext", realized_repository: "RealizedRepository", lint_ctx, lint_args):
    path = realized_repository.path
    for tool_path, tool_source in yield_tool_sources(ctx, path, recursive=True):
        original_path = tool_path.replace(path, realized_repository.real_path)
        info("+Linting tool %s" % original_path)
        if handle_tool_load_error(tool_path, tool_source):
            return True
        lint_tool_source_with(lint_ctx, tool_source, extra_modules=lint_args["extra_modules"])


def lint_shed_version(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.path

    tsi = tool_shed_instance("https://toolshed.g2.bx.psu.edu/")

    for tool_path, tool_source in yield_tool_sources(ctx=None, path=path, recursive=True):
        if handle_tool_load_error(tool_path, tool_source):
            continue

        if not isinstance(tool_source, ToolSource):
            continue

        repo_owner = realized_repository.owner
        repo_name = realized_repository.name
        tool_id = tool_source.parse_id()
        tool_version = parse_version(tool_source.parse_version())

        # check if there is already such a repo (otherwise get_ordered_installable_revisions will log an error for new repos)
        if len(tsi.repositories.get_repositories(repo_name, repo_owner)) == 0:
            continue

        try:
            installable_revisions = tsi.repositories.get_ordered_installable_revisions(repo_name, repo_owner)
        except ConnectionError:
            continue

        if len(installable_revisions) == 0:
            continue

        latest_installable_revision = installable_revisions[-1]
        repo_info, repo_metadata, _ = tsi.repositories.get_repository_revision_install_info(
            repo_name, repo_owner, latest_installable_revision
        )

        # no such tool in the TS -> fine to push any version
        if len(repo_metadata["valid_tools"]) == 0:
            continue

        # case 1 tool per repo
        if len(repo_metadata["valid_tools"]) == 1:
            assert repo_metadata["valid_tools"][0]["version"]
            ts_tool_version = repo_metadata["valid_tools"][0]["version"]
        # case n tools per repo
        else:
            tool = [_ for _ in repo_metadata["valid_tools"] if _["id"] == tool_id]
            assert len(tool) == 1
            assert tool[0]["version"]
            ts_tool_version = tool[0]["version"]

        if tool_version <= parse_version(ts_tool_version):
            lint_ctx.error(
                f"{tool_id}: version {tool_version} is less or equal than version of the latest installable revision {ts_tool_version}",
                "ShedVersion",
            )


def lint_expansion(realized_repository: "RealizedRepository", lint_ctx):
    missing = realized_repository.missing
    if missing:
        msg = "Failed to expand inclusions %s" % missing
        lint_ctx.warn(msg)
    else:
        lint_ctx.info("Included files all found.")


def lint_shed_metadata(realized_repository: "RealizedRepository", lint_ctx):
    found_all = True
    for key in SHED_METADATA:
        if key not in realized_repository.config:
            found_all = False
            lint_ctx.warn("Missing shed metadata field [%s] for repository" % key)
    if found_all:
        lint_ctx.info("Found all shed metadata fields required for automated repository creation and/or updates.")


def lint_readme(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.real_path
    readme_rst = os.path.join(path, "README.rst")
    readme = os.path.join(path, "README")
    readme_txt = os.path.join(path, "README.txt")

    readme_found = ""
    for readme in [readme_rst, readme, readme_txt]:
        if os.path.exists(readme):
            readme_found = readme

    if not readme_found:
        # TODO: filter on TYPE and make this a warning if
        # unrestricted repository - need to update iuc standards
        # first though.
        readme_md = os.path.join(path, "README.md")
        if os.path.exists(readme_md):
            lint_ctx.info("Found README in Markdown format, which is not rendered by the Tool Shed, skipping")
        else:
            lint_ctx.info("No README found, skipping.")
        return

    if readme_found.endswith(".rst"):
        with open(readme_found) as fh:
            readme_text = fh.read()
        invalid_rst = rst_invalid(readme_text)
        if invalid_rst:
            template = "Invalid restructured text found in README [%s]."
            msg = template % invalid_rst
            lint_ctx.warn(msg)
            return
        lint_ctx.info("README found containing valid reStructuredText.")
    else:
        lint_ctx.info("README found containing plain text.")


def lint_tool_dependencies_urls(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.real_path
    tool_dependencies = os.path.join(path, "tool_dependencies.xml")
    if not os.path.exists(tool_dependencies):
        lint_ctx.info("No tool_dependencies.xml, skipping.")
        return

    root = ET.parse(tool_dependencies).getroot()
    lint_urls(root, lint_ctx)


def lint_tool_dependencies_sha256sum(realized_repository: "RealizedRepository", lint_ctx):
    tool_dependencies = os.path.join(realized_repository.real_path, "tool_dependencies.xml")
    if not os.path.exists(tool_dependencies):
        lint_ctx.info("No tool_dependencies.xml, skipping.")
        return

    root = ET.parse(tool_dependencies).getroot()

    count = 0
    for action in root.findall(".//action"):
        assert action.tag == "action"
        if action.attrib.get("type", "") not in ["download_by_url", "download_file"]:
            continue
        url = (action.text or "").strip()
        checksum = action.attrib.get("sha256sum", "")
        if not checksum:
            lint_ctx.warn("Missing checksum for %s" % url)
        elif len(checksum) != 64 or not set("0123456789abcdef").issuperset(checksum.lower()):
            lint_ctx.error(f"Invalid checksum {checksum!r} for {url}")
        else:
            # TODO - See planned --verify option to check it matches
            # lint_ctx.info("SHA256 checkum listed for %s" % url)
            count += 1
    if count:
        lint_ctx.info("Found %i download action(s) with SHA256 checksums" % count)


def lint_tool_dependencies_xsd(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.real_path
    tool_dependencies = os.path.join(path, "tool_dependencies.xml")
    if not os.path.exists(tool_dependencies):
        lint_ctx.info("No tool_dependencies.xml, skipping.")
        return
    lint_xsd(lint_ctx, TOOL_DEPENDENCIES_XSD, tool_dependencies)


def lint_tool_dependencies_actions(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.real_path
    tool_dependencies = os.path.join(path, "tool_dependencies.xml")
    if not os.path.exists(tool_dependencies):
        lint_ctx.info("No tool_dependencies.xml, skipping.")
        return
    try:
        base.Dependencies(tool_dependencies)
        lint_ctx.info("Parsed tool dependencies.")
    except Exception as e:
        import sys
        import traceback

        exc_type, exc_value, exc_traceback = sys.exc_info()
        traceback.print_tb(exc_traceback, limit=1, file=sys.stdout)
        traceback.print_exc()
        template = "Problem parsing tool_dependenies.xml [%s]"
        msg = template % unicodify(e)
        lint_ctx.warn(msg)
        return


def lint_expected_files(realized_repository: "RealizedRepository", lint_ctx):
    if realized_repository.is_package:
        if not os.path.exists(realized_repository.tool_dependencies_path):
            lint_ctx.warn("Package repository does not contain a tool_dependencies.xml file.")

    if realized_repository.is_suite:
        if not os.path.exists(realized_repository.repo_dependencies_path):
            lint_ctx.warn("Suite repository does not contain a repository_dependencies.xml file.")


def lint_repository_dependencies(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.real_path
    repo_dependencies = os.path.join(path, "repository_dependencies.xml")
    if not os.path.exists(repo_dependencies):
        lint_ctx.info("No repository_dependencies.xml, skipping.")
        return
    lint_xsd(lint_ctx, REPO_DEPENDENCIES_XSD, repo_dependencies)


def lint_shed_yaml(realized_repository: "RealizedRepository", lint_ctx):
    path = realized_repository.real_path
    shed_yaml = os.path.join(path, ".shed.yml")
    if not os.path.exists(shed_yaml):
        lint_ctx.info("No .shed.yml file found, skipping.")
        return
    try:
        with open(shed_yaml) as fh:
            yaml.safe_load(fh)
    except Exception as e:
        lint_ctx.warn("Failed to parse .shed.yml file [%s]" % unicodify(e))
        return
    lint_ctx.info(".shed.yml found and appears to be valid YAML.")
    _lint_shed_contents(lint_ctx, realized_repository)


def _lint_shed_contents(lint_ctx, realized_repository: "RealizedRepository"):
    config = realized_repository.config

    def _lint_if_present(key, func, *args):
        value = config.get(key, None)
        if value is not None:
            msg = func(value, *args)
            if msg:
                lint_ctx.warn(msg)

    _lint_if_present("owner", validate_repo_owner)
    _lint_if_present("name", validate_repo_name)
    _lint_if_present("type", _validate_repo_type, config["name"])
    _lint_if_present("categories", _validate_categories, realized_repository)


def _validate_repo_type(repo_type, name):
    if repo_type not in VALID_REPOSITORY_TYPES:
        return "Invalid repository type specified [%s]" % repo_type

    is_dep = repo_type == "tool_dependency_definition"
    is_suite = repo_type == "repository_suite_definition"
    if is_dep and not name.startswith("package_"):
        return "Tool dependency definition repositories should have names starting with package_"
    if is_suite and not name.startswith("suite_"):
        return "Repository suite definition repositories should have names starting with suite_"
    if name.startswith("package_") or name.startswith("suite_"):
        if repo_type == "unrestricted":
            return "Repository name indicated specialized repository type but repository is listed as unrestricted."


def _validate_categories(categories, realized_repository: "RealizedRepository"):
    msg = None
    if len(categories) == 0:
        msg = "Repository should specify one or more categories."
    else:
        for category in categories:
            unknown_categories = []
            if category not in CURRENT_CATEGORIES:
                unknown_categories.append(category)
            if unknown_categories:
                msg = "Categories [%s] unknown." % unknown_categories
        if realized_repository.is_package:
            if "Tool Dependency Packages" not in categories:
                msg = "Packages should be placed and should only be placed in the category 'Tool Dependency Packages'."

    return msg


__all__ = ("lint_repository",)
