"""Module describing the planemo ``merge_test_reports`` command."""

import os

import click

from planemo import (
    io,
    options,
)
from planemo.cli import command_function
from planemo.galaxy.test.actions import merge_reports


@click.command("merge_test_reports")
@options.merge_test_json()
@options.tool_test_json("output_path")
@command_function
def cli(ctx, input_paths, output_path, **kwds):
    """Merge tool_test_output.json files from multiple runs."""
    for input_path in input_paths:
        if not os.path.exists(input_path):
            io.error("Failed to tool test json file at %s" % input_path)
            return 1
    merge_reports(input_paths=input_paths, output_path=output_path)
