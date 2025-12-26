# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json
import os
import sys

from spython.logger import bot
from spython.main.parse.parsers import get_parser
from spython.main.parse.writers import get_writer
from spython.utils import write_file, write_json


def main(args, options, parser):
    """This function serves as a wrapper around the DockerParser,
    SingularityParser, DockerWriter, and SingularityParser converters.
    We can either save to file if args.outfile is defined, or print
    to the console if not.
    """
    # We need something to work with
    if not args.files:
        parser.print_help()
        sys.exit(1)

    # Get the user specified input and output files
    outfile = None
    if len(args.files) > 1:
        outfile = args.files[1]

    # First try to get writer and parser, if not defined will return None
    writer = get_writer(args.writer)
    parser = get_parser(args.parser)

    # If the user wants to auto-detect the type
    if args.parser == "auto":
        if "dockerfile" in args.files[0].lower():
            parser = get_parser("docker")
        elif "singularity" in args.files[0].lower():
            parser = get_parser("singularity")

    # If the parser still isn't defined, no go.
    if parser is None:
        bot.exit(
            "Please provide a Dockerfile or Singularity recipe, or define the --parser type."
        )

    # If the writer needs auto-detect
    if args.writer == "auto":
        if parser.name == "docker":
            writer = get_writer("singularity")
        else:
            writer = get_writer("docker")

    # If the writer still isn't defined, no go
    if writer is None:
        bot.exit("Please define the --writer type.")

    # Initialize the chosen parser
    recipeParser = parser(args.files[0])

    # By default, discover entrypoint / cmd from Dockerfile
    entrypoint = "/bin/bash"
    force = False

    if args.entrypoint is not None:
        entrypoint = args.entrypoint

        # This is only done if the user intended to print json here
        recipeParser.entrypoint = args.entrypoint
        recipeParser.cmd = None
        force = True

    if args.json:
        if outfile is not None:
            if not os.path.exists(outfile):
                if force:
                    write_json(outfile, recipeParser.recipe.json())
                else:
                    bot.exit("%s exists, set --force to overwrite." % outfile)
        else:
            print(json.dumps(recipeParser.recipe.json(), indent=4))

    else:
        # Do the conversion
        recipeWriter = writer(recipeParser.recipe)
        result = recipeWriter.convert(runscript=entrypoint, force=force)

        # If the user specifies an output file, save to it
        if outfile is not None:
            write_file(outfile, result)

        # Otherwise, convert and print to screen
        else:
            print(result)
