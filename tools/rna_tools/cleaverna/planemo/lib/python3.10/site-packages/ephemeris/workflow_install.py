#!/usr/bin/env python
"""Tool to install workflows on a Galaxy instance."""
import argparse
import json
import os

from . import get_galaxy_connection
from .common_parser import (
    get_common_args,
    HideUnderscoresHelpFormatter,
)


def import_workflow(gi, path, publish_wf=False):
    """
    Given a connection to a Galaxy Instance (gi) and a path to a Galaxy workflow file,
    this function will import the worklfow into Galaxy.
    """
    with open(path) as wf_file:
        import_uuid = json.load(wf_file).get("uuid")
    existing_uuids = [d.get("latest_workflow_uuid") for d in gi.workflows.get_workflows()]
    if import_uuid not in existing_uuids:
        gi.workflows.import_workflow_from_local_path(path, publish=publish_wf)


def _parser():
    parent = get_common_args()
    parser = argparse.ArgumentParser(parents=[parent], formatter_class=HideUnderscoresHelpFormatter)
    parser.add_argument(
        "-w",
        "--workflow-path",
        "--workflow_path",
        required=True,
        help='Path to a workflow file or a directory with multiple workflow files ending with ".ga"',
    )
    parser.add_argument(
        "--publish-workflows",
        "--publish_workflows",
        required=False,
        action="store_true",
        help="Flag to publish all imported workflows, so that they are viewable by other users",
    )
    return parser


def main(argv=None):
    """
    This script uses bioblend to import .ga workflow files into a running instance of Galaxy
    """
    args = _parser().parse_args(argv)
    gi = get_galaxy_connection(args)

    if os.path.isdir(args.workflow_path):
        for file_path in os.listdir(args.workflow_path):
            if file_path.endswith(".ga"):
                import_workflow(
                    gi,
                    os.path.join(args.workflow_path, file_path),
                    args.publish_workflows,
                )
    else:
        import_workflow(gi, args.workflow_path, args.publish_workflows)


if __name__ == "__main__":
    main()
