#!/usr/bin/env python
"""Tool to generate tools from workflows"""
import json
from argparse import ArgumentParser
from collections.abc import Iterable

import yaml

from .common_parser import RawDescriptionHideUnderscoresHelpFormatter
from .shed_tools import InstallRepoDict
from .shed_tools_methods import format_tool_shed_url

INSTALL_TOOL_DEPENDENCIES = "install_tool_dependencies: True"
INSTALL_REPOSITORY_DEPENDENCIES = "install_repository_dependencies: True"
INSTALL_RESOLVER_DEPENDENCIES = "install_resolver_dependencies: True"


def _parser():
    parser = ArgumentParser(
        formatter_class=RawDescriptionHideUnderscoresHelpFormatter,
        usage="%(prog)s <options>",
        epilog="Workflow files must have been exported from Galaxy release 16.04 or newer.\n\n"
        "example:\n"
        "python %(prog)s -w workflow1 workflow2 -o mytool_list.yml -l my_panel_label\n"
        "Christophe Antoniewski <drosofff@gmail.com>\n"
        "https://github.com/ARTbio/ansible-artimed/tree/master/extra-files/generate_tool_list_from_ga_workflow_files.py",
    )
    parser.add_argument(
        "-w",
        "--workflow",
        dest="workflow_files",
        required=True,
        nargs="+",
        help="A space separated list of galaxy workflow description files in json format",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        required=True,
        dest="output_file",
        help="The output file with a yml tool list",
    )
    parser.add_argument(
        "-l",
        "--panel-label",
        "--panel_label",
        dest="panel_label",
        default="Tools from workflows",
        help="The name of the panel where the tools will show up in Galaxy." 'If not specified: "Tools from workflows"',
    )
    return parser


def get_workflow_dictionary(json_file):
    with open(json_file) as File:
        mydict = json.load(File)
    return mydict


def translate_workflow_dictionary_to_tool_list(workflow_dictionary, panel_label: str) -> list[InstallRepoDict]:
    starting_tool_list = extract_tool_shed_repositories_from_workflow_dict(workflow_dictionary)
    tool_list: list[InstallRepoDict] = []
    for tool in starting_tool_list:
        sub_dic: InstallRepoDict = {
            "name": tool["name"],
            "owner": tool["owner"],
            "revisions": [tool["changeset_revision"]],
            "tool_panel_section_label": panel_label,
            "tool_shed_url": format_tool_shed_url(tool["tool_shed"]),
        }
        tool_list.append(sub_dic)
    return tool_list


def extract_tool_shed_repositories_from_workflow_dict(workflow_dictionary):
    tool_list = []
    for step in workflow_dictionary["steps"].values():
        subworkflow = step.get("subworkflow")
        if subworkflow:
            tool_list.extend(extract_tool_shed_repositories_from_workflow_dict(subworkflow))
        tsr = step.get("tool_shed_repository")
        if tsr:
            tool_list.append(tsr)
    return tool_list


def print_yaml_tool_list(tool_dictionary, output_file):
    with open(output_file, "w") as F:
        F.write(
            "\n".join(
                [
                    INSTALL_TOOL_DEPENDENCIES,
                    INSTALL_REPOSITORY_DEPENDENCIES,
                    INSTALL_RESOLVER_DEPENDENCIES,
                    "",
                    "",
                ]
            )
        )
        F.write(yaml.safe_dump(tool_dictionary, default_flow_style=False))
    return


def reduce_tool_list(tool_list: list[InstallRepoDict]) -> list[InstallRepoDict]:
    for current_tool in tool_list:
        for tool in tool_list:
            if current_tool is tool:
                continue
            if (
                tool["name"] == current_tool["name"]
                and tool["owner"] == current_tool["owner"]
                and tool["tool_panel_section_label"] == current_tool["tool_panel_section_label"]
                and tool["tool_shed_url"] == current_tool["tool_shed_url"]
            ):
                current_tool["revisions"].extend(tool["revisions"])
                tool_list.remove(tool)
        current_tool["revisions"] = list(set(current_tool["revisions"]))
    return tool_list


def generate_repo_list_from_workflow(workflow_files: Iterable[str], panel_label: str) -> list[InstallRepoDict]:
    intermediate_tool_list: list[InstallRepoDict] = []
    for workflow in workflow_files:
        workflow_dictionary = get_workflow_dictionary(workflow)
        intermediate_tool_list += translate_workflow_dictionary_to_tool_list(workflow_dictionary, panel_label)
    return reduce_tool_list(intermediate_tool_list)


def generate_tool_list_from_workflow(workflow_files: Iterable[str], panel_label: str, output_file: str):
    """
    :rtype: object
    """

    convert_dict = {"tools": generate_repo_list_from_workflow(workflow_files=workflow_files, panel_label=panel_label)}
    print_yaml_tool_list(convert_dict, output_file)


def main(argv=None):
    options = _parser().parse_args(argv)
    generate_tool_list_from_workflow(options.workflow_files, options.panel_label, options.output_file)


if __name__ == "__main__":
    main()
