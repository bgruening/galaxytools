#!/usr/bin/env python
"""Tool to install tool dependencies on a Galaxy instance."""
import argparse
import logging as log
import os
import xml.etree.ElementTree as ET

import yaml
from bioblend import ConnectionError as ConnErr
from bioblend.galaxy.tools import ToolClient

from ephemeris import get_galaxy_connection
from ephemeris.common_parser import (
    get_common_args,
    HideUnderscoresHelpFormatter,
)

timeout_codes = (408, 502, 504)


def _parser():
    parent = get_common_args()
    parser = argparse.ArgumentParser(parents=[parent], formatter_class=HideUnderscoresHelpFormatter)
    parser.add_argument(
        "-t",
        "--tool",
        help="Path to a tool file, tool_conf file, or yaml file containing a sequence of tool ids",
        nargs="*",
    )
    parser.add_argument("-i", "--id", help="Space-separated list of tool ids", nargs="*")

    return parser


def _install(tool_client, tool_id):
    try:
        tool_client.install_dependencies(tool_id)
    except ConnErr as e:
        if e.status_code in timeout_codes:
            log.warning(e.body)
        else:
            raise


def main(argv=None):
    """
    This script uses bioblend to trigger dependencies installations for the provided tools
    """
    args = _parser().parse_args(argv)
    gi = get_galaxy_connection(args)
    tool_client = ToolClient(gi)

    if args.verbose:
        log.basicConfig(level=log.DEBUG)

    if args.tool:
        for tool_conf_path in args.tool:  # type: str
            _, ext = os.path.splitext(tool_conf_path)
            if ext == ".xml":
                log.info("tool_conf xml found, parsing..")
                # install all
                root = ET.ElementTree(file=tool_conf_path).getroot()
                if root.tag == "toolbox":
                    # Install all from tool_conf
                    tool_path = root.get("tool_path", "")
                    tool_path = tool_path.replace(
                        "${tool_conf_dir}",
                        os.path.abspath(os.path.dirname(tool_conf_path)),
                    )
                    if tool_path:
                        log.info("Searching for tools relative to %s", tool_path)
                    tools = root.findall(".//tool[@file]")
                    if len(tools) == 0:
                        log.warning("No tools found in tool_conf")
                        continue

                    for tool in tools:
                        tool_id = ET.ElementTree(file=os.path.join(tool_path, tool.get("file"))).getroot().get("id")
                        if tool_id:
                            log.info(
                                "Installing tool dependencies for %s from: %s",
                                tool_id,
                                tool.get("file"),
                            )
                            _install(tool_client, tool_id)
                elif root.tag == "tool" and root.get("id"):
                    # Install from single tool file
                    log.info("Tool xml found. Installing %s dependencies", root.get("id"))
                    _install(tool_client, root.get("id"))
            else:
                log.info("YAML tool list found, parsing..")
                with open(tool_conf_path) as fh:
                    tool_ids = yaml.safe_load(fh)
                for tool_id in tool_ids:
                    # Install from yaml file
                    log.info("Installing %s dependencies..", tool_id)
                    _install(tool_client, tool_id)

    if args.id:
        for tool_id in args.id:  # type: str
            log.info("Installing %s dependencies..", tool_id)
            _install(tool_client, tool_id.strip())


if __name__ == "__main__":
    main()
