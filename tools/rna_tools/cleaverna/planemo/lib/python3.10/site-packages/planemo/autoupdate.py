"""Autoupdate older conda dependencies in the requirements section."""

import collections
import itertools
import re
import xml.etree.ElementTree as ET
from typing import (
    Any,
    DefaultDict,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    TYPE_CHECKING,
    Union,
)
from xml.etree.ElementTree import (
    Element,
    ElementTree,
)

import requests
import yaml
from bioblend.toolshed import ToolShedInstance
from galaxy.tool_util.deps import conda_util
from galaxy.tool_util.version import parse_version
from packaging.version import Version

import planemo.conda
from planemo.galaxy.workflows import (
    get_tool_ids_for_workflow,
    get_toolshed_url_for_tool_id,
)
from planemo.io import (
    error,
    info,
)
from planemo.workflow_lint import find_repos_from_tool_id

if TYPE_CHECKING:
    from planemo.cli import PlanemoCliContext
    from planemo.galaxy.config import LocalGalaxyConfig
    from planemo.runnable import Runnable


def bump_version(version_str: str) -> str:
    version = Version(version_str)

    # Extract base version and pre-release parts
    base_version = list(version.release)
    pre = version.pre
    post = version.post
    dev = version.dev

    if pre:
        # Increment the pre-release number (e.g., beta1 -> beta2)
        new_pre = (pre[0], pre[1] + 1)
        return f"{'.'.join(map(str, base_version))}{new_pre[0]}{new_pre[1]}"
    elif post:
        return f"{version.base_version}.post{post + 1}"
    elif dev:
        return f"{version.base_version}.dev{dev + 1}"
    else:
        # Find the most minor specified version part and increment it
        base_version[-1] += 1
        return ".".join(map(str, base_version))


def find_macros(xml_tree: "ElementTree[Element[str]]") -> List[Any]:
    """
    Get macros from the XML tree
    """
    macros = []
    for macro_import in xml_tree.iter("import"):
        macros.append(macro_import.text)
    return macros


def get_requirements(
    xml_tree: "ElementTree[Element[str]]",
) -> Tuple[Dict[str, Dict[str, Optional[str]]], Optional[str]]:
    """
    Get requirements from the XML tree
    """
    requirements = {}
    main_req = None
    for requirement in xml_tree.iter("requirement"):
        if requirement.attrib.get("version") == "@TOOL_VERSION@":
            main_req = requirement.text
        else:
            assert requirement.text
            requirements[requirement.text] = {
                "tag": ET.tostring(requirement, encoding="unicode").strip(),
                "text": requirement.attrib.get("version"),
            }
    return requirements, main_req


def get_tokens(xml_tree: "ElementTree[Element[str]]") -> Dict[str, Dict[str, Optional[str]]]:
    """
    Get tokens from the XML tree
    """
    tokens = {}
    for token in xml_tree.iter("token"):
        tokens[token.attrib["name"]] = {"tag": ET.tostring(token, encoding="unicode").strip(), "text": token.text}
    return tokens


def check_conda(package_name: str, ctx: "PlanemoCliContext", **kwds) -> str:
    """
    Get the most up-to-date conda version for a package.
    """
    conda_context = planemo.conda.build_conda_context(ctx, **kwds)
    if not conda_context.is_conda_installed():
        # check directly via Anaconda API
        r = requests.get("https://api.anaconda.org/search", params={"name": package_name})
        search_results = itertools.chain.from_iterable(
            n["versions"] for n in r.json() if n["name"] == package_name and n["owner"] in kwds["conda_ensure_channels"]
        )
        return sorted(search_results, key=parse_version, reverse=True)[0]

    target = conda_util.CondaTarget(package_name)
    best_search_results = conda_util.best_search_result(target, conda_context=conda_context)
    if best_search_results[0] is None:
        raise Exception(f"No conda package found for {package_name}")
    return best_search_results[0]["version"]


def update_xml(
    tool_path: str,
    xml_tree: "ElementTree[Element[str]]",
    tags_to_update: List[Dict[str, str]],
    wrapper_version_token: Optional[Union[int, str]],
    is_macro: bool = False,
) -> None:
    """
    Write modified XML to tool_path
    """

    def update_token(xml_text, tag, token_value):
        new_tag = f">{token_value}<".join(re.split(">.*<", tag))
        return re.sub(tag, new_tag, xml_text)

    def update_requirement(xml_text, tag, requirement_value):
        new_tag = f'version="{requirement_value}"'.join(re.split('version=".*"', tag))
        return re.sub(tag, new_tag, xml_text)

    with open(tool_path, "r+", newline="") as f:
        xml_text = f.read()
        for tag_to_update in tags_to_update:
            if tag_to_update["type"] == "token":
                xml_text = update_token(xml_text, tag_to_update["tag"], tag_to_update["value"])
            if tag_to_update["type"] == "requirement":
                xml_text = update_requirement(xml_text, tag_to_update["tag"], tag_to_update["value"])
        if wrapper_version_token == 0 and not is_macro:
            # i.e. @VERSION_SUFFIX@ not specified so update the version directly in the tool tag
            tool_tag = re.sub(
                'version="@TOOL_VERSION@.*?"',
                'version="@TOOL_VERSION@+galaxy0"',
                re.findall('<tool .*version="@TOOL_VERSION@.*">', xml_text)[0],
            )
            xml_text = re.sub('<tool .*version="@TOOL_VERSION@.*">', tool_tag, xml_text)
        f.seek(0)
        f.truncate()
        f.write(xml_text)


def create_requirement_dict(
    xml_files: Dict[str, "ElementTree[Element[str]]"], skip_reqs: List[str]
) -> Tuple[Dict[str, Dict[str, Dict[str, Optional[str]]]], Optional[Tuple[str, str]]]:
    """
    Create dict with requirements and find main requirement
    """
    requirements = {}
    main_req = None
    for k, v in xml_files.items():
        file_reqs, file_main_req = get_requirements(v)
        requirements[k] = {k: v for k, v in file_reqs.items() if k not in skip_reqs}
        if file_main_req:
            if main_req:
                error("Multiple requirements use the token @TOOL_VERSION@!")
            main_req = (file_main_req, k)
    if not main_req:
        error("No requirement uses the token @TOOL_VERSION@!")
    return requirements, main_req


def create_token_dict(
    ctx: "PlanemoCliContext", xml_files: Dict[str, "ElementTree[Element[str]]"], main_req: Tuple[str, str], **kwds
) -> Tuple[
    Dict[str, Dict[str, Dict[str, Optional[str]]]], DefaultDict[str, List[Dict[str, str]]], Optional[str], Optional[str]
]:
    """
    Create dict with relevant tokens and check conda requirements for main
    """
    tokens: Dict[str, Dict[str, Dict[str, Optional[str]]]] = {}
    current_main_req, updated_main_req = None, None
    xml_to_update = collections.defaultdict(list)
    for k, v in xml_files.items():
        tokens[k] = get_tokens(v)
        # check if it is @TOOL_VERSION@ and if so do check_conda
        if "@TOOL_VERSION@" in tokens[k]:
            current_main_req = tokens[k]["@TOOL_VERSION@"]["text"]
            updated_main_req = check_conda(main_req[0], ctx, **kwds)
            if current_main_req:
                tag = tokens[k]["@TOOL_VERSION@"]["tag"]
                assert tag is not None
                xml_to_update[k].append({"type": "token", "tag": tag, "value": updated_main_req})

    return tokens, xml_to_update, current_main_req, updated_main_req


def perform_required_update(
    ctx: "PlanemoCliContext",
    xml_files: Dict[str, "ElementTree[Element[str]]"],
    tool_path: str,
    requirements: Dict[str, Dict[str, Dict[str, Optional[str]]]],
    tokens: Dict[str, Dict[str, Dict[str, Optional[str]]]],
    xml_to_update: DefaultDict[str, List[Dict[str, str]]],
    wrapper_version_token: Optional[Union[int, str]],
    **kwds,
) -> Set[str]:
    """
    Carry out the update, if requirements are out-of-date
    """
    # check all requirements
    for k, v in requirements.items():
        for req in v:
            req_check = check_conda(req, ctx, **kwds)
            # print(req_check, v[req]['text'])
            if req_check != v[req]["text"]:
                tag = v[req]["tag"]
                assert tag is not None
                xml_to_update[k].append({"type": "requirement", "tag": tag, "value": req_check})

    # check all tokens, if wrapper_version_token exists
    if wrapper_version_token:
        for k, v in tokens.items():
            if isinstance(wrapper_version_token, str) and wrapper_version_token in v:
                tag = v[wrapper_version_token]["tag"]
                assert tag is not None
                xml_to_update[k].append({"type": "token", "tag": tag, "value": "0"})

    # finally, update each file separately
    for k, et in xml_files.items():
        update_xml(k, et, xml_to_update[k], wrapper_version_token, is_macro=(k != tool_path))
    info(f"Tool {tool_path} successfully updated.")
    return set(xml_files)


def autoupdate_tool(ctx: "PlanemoCliContext", tool_path: str, modified_files: Set[Any], **kwds) -> Optional[Set[str]]:
    """
    Autoupdate an XML file
    """
    modified_files = modified_files or set()
    # create a dict of all files that need editing - wrapper plus macros
    xml_files = {tool_path: ET.parse(tool_path)}

    # get name of token which defines the wrapper version; if just an integer, None
    version_str = xml_files[tool_path].getroot().attrib.get("version")
    if version_str:
        versions = version_str.split("+galaxy")
        if versions[0] != "@TOOL_VERSION@":
            error("Tool version does not contain @TOOL_VERSION@ as required by autoupdate.")
            return None
        elif len(versions) == 1:
            wrapper_version_token = None
        else:
            if versions[1][0] == versions[1][-1] == "@":
                wrapper_version_token = versions[1]
            else:
                wrapper_version_token = "0"  # assume an int, reset to 0
    else:
        wrapper_version_token = None

    # add macros to xml_files
    macro_paths = xml_files[tool_path]
    if macro_paths:
        for macro in find_macros(macro_paths):
            macro_path = "/".join(tool_path.split("/")[:-1] + [macro])
            xml_files[macro_path] = ET.parse(macro_path)

    requirements, main_req = create_requirement_dict(xml_files, kwds.get("skip_requirements", "").split(","))
    if main_req is None:
        return None
    tokens, xml_to_update, current_main_req, updated_main_req = create_token_dict(ctx, xml_files, main_req, **kwds)

    if current_main_req == updated_main_req and not (modified_files & set(xml_files)):
        info(f"No updates required or made to {tool_path}.")
        return None  # end here if no update needed

    if kwds.get("dry_run"):
        error(
            f"Update required to {tool_path}! Tool main requirement has version {current_main_req}, newest conda version is {updated_main_req}"
        )
        return None

    else:
        info(f"Updating {tool_path.split('/')[-1]} from version {current_main_req} to {updated_main_req}")
        return perform_required_update(
            ctx, xml_files, tool_path, requirements, tokens, xml_to_update, wrapper_version_token, **kwds
        )


def _update_wf(config: "LocalGalaxyConfig", workflow_id: str, instance: bool = False) -> None:
    """
    Recursively update a workflow, including subworkflows
    """
    wf = config.user_gi.make_get_request(
        f"{config.user_gi.url}/workflows/{workflow_id}", params={"instance": instance}
    ).json()
    for step in wf.get("steps", {}).values():
        if step["type"] == "subworkflow":
            # update subworkflows before the main workflow
            _update_wf(config, step["workflow_id"], instance=True)
    config.user_gi.workflows.refactor_workflow(wf["id"], actions=[{"action_type": "upgrade_all_steps"}])


def get_newest_tool_id(tool_ids: List[str]) -> str:
    return sorted(
        tool_ids,
        key=lambda n: parse_version(n.split("/")[-1]),
    )[-1]


def outdated_tools(  # noqa: C901
    ctx: "PlanemoCliContext", wf_dict: Dict[str, Any], tools_to_skip: List[str]
) -> Dict[str, Dict[str, str]]:
    """
    tools_to_skip should be a list of base tool ids.
    """

    def base_tool_id(tool_id: str) -> str:
        return tool_id.rsplit("/", 1)[0]

    def check_tool_step(tool_id: str) -> Dict[str, Dict[str, str]]:
        """
        Return a dict with current and newest tool version, in case they don't match.
        """
        tool_shed_url = get_toolshed_url_for_tool_id(tool_id)
        if not tool_shed_url:
            return {}
        ts = ToolShedInstance(tool_shed_url)
        warning_msg, repos = find_repos_from_tool_id(tool_id, ts)
        if warning_msg != "":
            ctx.log(warning_msg)
        if len(repos) == 0:
            return repos
        base_id = base_tool_id(tool_id)
        matching_tool_ids = []
        for repo in repos.values():
            if isinstance(repo, dict):
                for tool in repo.get("tools") or []:
                    if tool["guid"].startswith(base_id + "/"):
                        matching_tool_ids.append(tool["guid"])
                        # there can only be one matching tool id in a repo
                        break
        updated_tool_id = get_newest_tool_id(matching_tool_ids)
        if tool_id != updated_tool_id:
            return {base_id: {"current": tool_id, "updated": updated_tool_id}}
        else:
            return {}

    def outdated_tools_rec(wf_dict: Dict[str, Any]) -> None:
        tool_ids = get_tool_ids_for_workflow(wf_dict)
        for tool_id in tool_ids:
            base_id = base_tool_id(tool_id)
            if base_id not in checked_tools:
                outdated_tool_dict.update(check_tool_step(tool_id))
                checked_tools.append(base_id)

    outdated_tool_dict: Dict[str, Dict[str, str]] = {}
    # Initialize the list of tools already checked with a copy of tools_to_skip
    checked_tools = tools_to_skip.copy()
    outdated_tools_rec(wf_dict)
    return outdated_tool_dict


def get_tools_to_update(
    ctx: "PlanemoCliContext", workflow: "Runnable", tools_to_skip: List[str]
) -> Dict[str, Dict[str, str]]:
    # before we run the autoupdate, we check the tools against the toolshed to see if there
    # are any new versions. This saves spinning up Galaxy and installing the tools if there
    # is nothing to do, and also allows us to collect a list of the tools which need updating
    with open(workflow.path) as f:
        wf_dict = yaml.load(f, Loader=yaml.SafeLoader)

    return outdated_tools(ctx, wf_dict, tools_to_skip)


def autoupdate_wf(ctx: "PlanemoCliContext", config: "LocalGalaxyConfig", wf: "Runnable") -> Dict[str, Any]:
    workflow_id = config.workflow_id_for_runnable(wf)
    _update_wf(config, workflow_id)
    return config.user_gi.workflows.export_workflow_dict(workflow_id)


def fix_workflow(original_wf: Dict[str, Any], updated_wf: Dict[str, Any]) -> Dict[str, Any]:
    # the Galaxy refactor action can't do everything right now... some manual fixes here
    # * bump release number if present
    # check release; bump if it exists
    if original_wf.get("release"):
        updated_wf["release"] = bump_version(original_wf["release"])
    return updated_wf
