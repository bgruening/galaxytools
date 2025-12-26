#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
import argparse
import os
import sys
from collections.abc import Iterator
from typing import cast

import ruamel.yaml

import cwl_utils.parser as cwl
from cwl_utils.image_puller import (
    DockerImagePuller,
    ImagePuller,
    SingularityImagePuller,
)


def arg_parser() -> argparse.ArgumentParser:
    """Argument parser."""
    parser = argparse.ArgumentParser(
        description="Save container images specified in a CWL document (Workflow or CommandLineTool). "
        "For CWL Workflows, all steps will also be searched (recursively)."
    )
    parser.add_argument(
        "input", help="Input CWL document (CWL Workflow or CWL CommandLineTool)"
    )
    parser.add_argument("--dir", help="Directory in which to save images")
    parser.add_argument(
        "-s",
        "--singularity",
        help="Use singularity to pull the image",
        action="store_true",
    )
    parser.add_argument(
        "--container-engine",
        dest="container_engine",
        help="Specify which command to use to run OCI containers. "
        "Defaults to 'docker' (or 'singularity' if --singularity/-s is passed).",
    )
    parser.add_argument(
        "--force-download", help="Force pulling a newer container.", action="store_true"
    )
    return parser


def run(args: argparse.Namespace) -> list[cwl.DockerRequirement]:
    """Extract the docker reqs and download them using Singularity or Docker."""
    if args.singularity and not args.dir:
        print("Error! Must specify --dir if using --singularity")
        sys.exit(1)

    if args.dir:
        os.makedirs(args.dir, exist_ok=True)

    top = cwl.load_document_by_uri(args.input)
    reqs: list[cwl.DockerRequirement] = []

    for req in traverse(top):
        reqs.append(req)
        if not req.dockerPull:
            print(
                "Unable to save image from due to lack of 'dockerPull':",
                file=sys.stderr,
            )
            yaml = ruamel.yaml.YAML()
            yaml.dump(req.save(), sys.stderr)
            continue
        if args.singularity:
            image_puller: ImagePuller = SingularityImagePuller(
                req.dockerPull,
                args.dir,
                (
                    args.container_engine
                    if args.container_engine is not None
                    else "singularity"
                ),
                args.force_download,
            )
        else:
            image_puller = DockerImagePuller(
                req.dockerPull,
                args.dir,
                (
                    args.container_engine
                    if args.container_engine is not None
                    else "docker"
                ),
                args.force_download,
            )
        image_puller.save_docker_image()
    return reqs


def extract_docker_requirements(
    process: cwl.Process,
) -> Iterator[cwl.DockerRequirement]:
    """Yield an iterator of the docker reqs, normalizing the pull request."""
    for req in extract_docker_reqs(process):
        if isinstance(req.dockerPull, str) and ":" not in req.dockerPull:
            req.dockerPull += ":latest"
        yield req


def extract_docker_reqs(process: cwl.Process) -> Iterator[cwl.DockerRequirement]:
    """For the given process, extract the DockerRequirement(s)."""
    if process.requirements:
        for req in process.requirements:
            if isinstance(req, cwl.DockerRequirementTypes):
                yield req
    if process.hints:
        for req in process.hints:
            if isinstance(req, cwl.DockerRequirementTypes):
                yield req


def traverse(process: cwl.Process) -> Iterator[cwl.DockerRequirement]:
    """Yield the iterator for the docker reqs, including an workflow steps."""
    yield from extract_docker_requirements(process)
    if isinstance(process, cwl.WorkflowTypes):
        yield from traverse_workflow(process)


def get_process_from_step(step: cwl.WorkflowStep) -> cwl.Process:
    """Return the process for this step, loading it if necessary."""
    if isinstance(step.run, str):
        return cast(cwl.Process, cwl.load_document_by_uri(step.run))
    return cast(cwl.Process, step.run)


def traverse_workflow(workflow: cwl.Workflow) -> Iterator[cwl.DockerRequirement]:
    """Iterate over the steps of this workflow, yielding the docker reqs."""
    for step in workflow.steps:
        yield from extract_docker_reqs(step)
        yield from traverse(get_process_from_step(step))


def main() -> int:
    """Command line entry point."""
    run(arg_parser().parse_args(sys.argv[1:]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
