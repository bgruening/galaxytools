#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright 2021 Michael R. Crusoe
"""Normalize CWL documents to CWL v1.2, JSON style."""
import argparse
import logging
import sys
import tempfile
from collections.abc import MutableSequence
from pathlib import Path

from cwlupgrader import main as cwlupgrader
from ruamel import yaml
from schema_salad.sourceline import add_lc_filename

from cwl_utils import cwl_v1_2_expression_refactor
from cwl_utils.loghandler import _logger as _cwlutilslogger
from cwl_utils.pack import pack
from cwl_utils.parser.cwl_v1_2 import load_document_by_yaml, save

_logger = logging.getLogger("cwl-normalizer")  # pylint: disable=invalid-name
defaultStreamHandler = logging.StreamHandler()  # pylint: disable=invalid-name
_logger.addHandler(defaultStreamHandler)
_logger.setLevel(logging.INFO)
_cwlutilslogger.setLevel(100)


def arg_parser() -> argparse.ArgumentParser:
    """Build the argument parser."""
    parser = argparse.ArgumentParser(
        description="Normalizes CWL documents. Will upgrade to CWL v1.2, "
        "and pack the result. Can optionally refactor out CWL expressions."
    )
    parser.add_argument(
        "--etools",
        help="Output ExpressionTools, don't go all the way to CommandLineTools.",
        action="store_true",
    )
    parser.add_argument(
        "--skip-some1",
        help="Don't process CommandLineTool.inputs.inputBinding and CommandLineTool.arguments sections.",
        action="store_true",
    )
    parser.add_argument(
        "--skip-some2",
        help="Don't process CommandLineTool.outputEval or "
        "CommandLineTool.requirements.InitialWorkDirRequirement.",
        action="store_true",
    )
    parser.add_argument(
        "--no-expression-refactoring",
        help="Don't do any CWL expression refactoring.",
        action="store_true",
    )
    parser.add_argument("dir", help="Directory in which to save converted files")
    parser.add_argument(
        "inputs",
        nargs="+",
        help="One or more CWL documents.",
    )
    return parser


def parse_args(args: list[str]) -> argparse.Namespace:
    """Parse the command line arguments."""
    return arg_parser().parse_args(args)


def main() -> None:
    """Console entry point."""
    sys.exit(run(parse_args(sys.argv[1:])))


def run(args: argparse.Namespace) -> int:
    """Primary processing loop."""
    imports: set[str] = set()
    for document in args.inputs:
        _logger.info("Processing %s.", document)
        with open(document) as doc_handle:
            result = yaml.main.round_trip_load(doc_handle, preserve_quotes=True)
        add_lc_filename(result, document)
        version = result.get("cwlVersion", None)
        if version in ("draft-3", "cwl:draft-3", "v1.0", "v1.1"):
            result = cwlupgrader.upgrade_document(result, args.dir, imports=imports)
        else:
            _logger.error(
                "Sorry, %s in %s is not a supported CWL version by this tool.",
                (version, document),
            )
            return -1
        uri = Path(document).resolve().as_uri()
        if not args.no_expression_refactoring:
            refactored, _ = cwl_v1_2_expression_refactor.traverse(
                load_document_by_yaml(result, uri),
                not args.etools,
                False,
                args.skip_some1,
                args.skip_some2,
            )
            if not isinstance(refactored, MutableSequence):
                result = save(
                    refactored,
                    base_url=(
                        refactored.loadingOptions.fileuri
                        if refactored.loadingOptions.fileuri
                        else ""
                    ),
                )
            #   ^^ Setting the base_url and keeping the default value
            #      for relative_uris=True means that the IDs in the generated
            #      JSON/YAML are kept clean of the path to the input document
            else:
                result = [
                    save(result_item, base_url=result_item.loadingOptions.fileuri)
                    for result_item in refactored
                ]
        if "$graph" in result:
            packed = result
        else:
            with tempfile.TemporaryDirectory() as tmpdirname:
                path = Path(tmpdirname) / Path(document).name
                packed = pack(str(path))
        output = Path(args.dir) / Path(document).name
        with open(output, "w", encoding="utf-8") as output_filehandle:
            output_filehandle.write(packed)
    return 0


if __name__ == "__main__":
    main()
