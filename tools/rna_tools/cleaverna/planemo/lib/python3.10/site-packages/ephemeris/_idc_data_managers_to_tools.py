#!/usr/bin/env python
"""Helper script for IDC - not yet meant for public consumption.

This script takes a data_managers.yml configuration describing the
set of data managers the IDC configuration targets and builds a
a tools.yml file from it for use with shed_tools.
"""
import argparse
import logging
from typing import NamedTuple

import yaml

from ._config_models import (
    read_data_managers,
    RepositoryInstallTargets,
)
from .common_parser import (
    add_log_file_argument,
    add_verbosity_argument,
)
from .ephemeris_log import (
    disable_external_library_logging,
    setup_global_logger,
)


class DataManager(NamedTuple):
    tool_id: str
    repository_name: str
    tags: list[str]


def read_data_managers_configuration(path: str) -> dict[str, DataManager]:
    raw_data_managers = read_data_managers(path)
    data_managers: dict[str, DataManager] = {}
    for repository_name, data_manager_configuration in raw_data_managers.root.items():
        data_manager = DataManager(
            tool_id=data_manager_configuration.tool_id,
            repository_name=repository_name,
            tags=data_manager_configuration.tags or [],
        )
        data_managers[repository_name] = data_manager
    return data_managers


def build_shed_install_conf(path: str) -> dict:
    data_managers = read_data_managers_configuration(path)
    tools = []
    for data_manager in data_managers.values():
        tool_id = data_manager.tool_id
        tool_id_parts = tool_id.split("/")
        repo_owner = tool_id_parts[2]
        repo_name = tool_id_parts[3]
        entry = {
            "name": repo_name,
            "owner": repo_owner,
            "tool_panel_section_label": None,
            "tool_shed_url": "toolshed.g2.bx.psu.edu",
        }
        tools.append(entry)
    tools_yaml = {"tools": tools}
    return tools_yaml


def write_shed_install_conf(data_manager_conf_path: str, output_path: str) -> None:
    tools_yaml = build_shed_install_conf(data_manager_conf_path)

    # validate generated dict to ensure we're writing out valid file
    RepositoryInstallTargets(**tools_yaml)

    with open(output_path, "w") as f:
        yaml.safe_dump(tools_yaml, f)


def _parser():
    """returns the parser object."""

    parser = argparse.ArgumentParser(add_help=False)
    general_group = parser.add_argument_group("General options")
    add_verbosity_argument(general_group)
    add_log_file_argument(general_group)
    parser.add_argument("--data-managers-conf", default="data_managers.yml")
    parser.add_argument("--shed-install-output-conf", default="tools.yml")
    return parser


def main():
    disable_external_library_logging()
    parser = _parser()
    args = parser.parse_args()
    log = setup_global_logger(name=__name__, log_file=args.log_file)
    if args.verbose:
        log.setLevel(logging.DEBUG)
    else:
        log.setLevel(logging.INFO)
    write_shed_install_conf(args.data_managers_conf, args.shed_install_output_conf)


if __name__ == "__main__":
    main()
