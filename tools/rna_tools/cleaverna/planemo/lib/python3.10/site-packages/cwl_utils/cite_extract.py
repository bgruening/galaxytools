#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import argparse
import sys
from collections.abc import Iterator
from typing import cast

import cwl_utils.parser as cwl


def arg_parser() -> argparse.ArgumentParser:
    """Construct the argument parser."""
    parser = argparse.ArgumentParser(
        description="Print information about software used in a CWL document (Workflow or CommandLineTool). "
        "For CWL Workflows, all steps will also be searched (recursively)."
    )
    parser.add_argument(
        "input", help="Input CWL document (CWL Workflow or CWL CommandLineTool)"
    )
    return parser


def run(args: argparse.Namespace) -> int:
    """Extract the software requirements."""
    for req in traverse(cwl.load_document_by_uri(args.input)):
        process_software_requirement(req)
    return 0


def main() -> int:
    """Console entry point."""
    return run(arg_parser().parse_args(sys.argv[1:]))


def extract_software_reqs(
    process: cwl.Process,
) -> Iterator[cwl.SoftwareRequirement]:
    """Return an iterator over any SoftwareRequirements found in the given process."""
    if process.requirements:
        for req in process.requirements:
            if isinstance(req, cwl.SoftwareRequirementTypes):
                yield req
    if process.hints:
        for req in process.hints:
            if isinstance(req, cwl.SoftwareRequirementTypes):
                yield req


def process_software_requirement(req: cwl.SoftwareRequirement) -> None:
    """Pretty print the software package information."""
    for package in req.packages:
        print(
            "Package: {}, version: {}, specs: {}".format(
                package.package, package.version, package.specs
            )
        )


def traverse(process: cwl.Process) -> Iterator[cwl.SoftwareRequirement]:
    """Extract the software packages for this process, and any steps."""
    yield from extract_software_reqs(process)
    if isinstance(process, cwl.WorkflowTypes):
        yield from traverse_workflow(process)


def get_process_from_step(step: cwl.WorkflowStep) -> cwl.Process:
    """Return the process for this step, loading it if needed."""
    if isinstance(step.run, str):
        return cast(cwl.Process, cwl.load_document_by_uri(step.run))
    return cast(cwl.Process, step.run)


def traverse_workflow(workflow: cwl.Workflow) -> Iterator[cwl.SoftwareRequirement]:
    """Iterate over the given workflow, extracting the software packages."""
    for step in workflow.steps:
        yield from extract_software_reqs(step)
        yield from traverse(get_process_from_step(step))


if __name__ == "__main__":
    sys.exit(main())
