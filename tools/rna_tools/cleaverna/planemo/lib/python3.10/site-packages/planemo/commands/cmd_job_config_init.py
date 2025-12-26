"""Module describing the planemo ``project_init`` command."""

import os
import sys

import click
from gxjobconfinit import (
    build_job_config,
    ConfigArgs,
)
from gxjobconfinit.generate import DevelopmentContext

from planemo import options
from planemo.cli import command_function
from planemo.galaxy.config import get_all_tool_path_from_kwds
from planemo.io import (
    info,
    launch_if_open_flagged,
    warn,
)
from planemo.runnable import for_paths
from planemo.tools import uris_to_paths

SUCCESS_MESSAGE = "Wrote configuration template to %s, please open with editor to customize if needed."


@click.command("job_config_init")
@options.optional_tools_arg(multiple=True, allow_uris=True)
@options.open_file_option()
@options.job_config_init_options()
@command_function
def cli(ctx, uris, **kwds):
    """Initialize an small Galaxy job config file for hosted workflow runs."""
    config_path = "job_conf.yml"
    if os.path.exists(config_path):
        warn(f"File '{config_path}' already exists, exiting.")
        sys.exit(1)

    paths = uris_to_paths(ctx, uris)
    runnables = for_paths(paths)
    tool_paths = get_all_tool_path_from_kwds(runnables, **kwds)
    dev_context = DevelopmentContext(
        kwds.get("test_data", None),
        tool_paths,
    )
    init_config = ConfigArgs.from_dict(**kwds)
    job_config = build_job_config(init_config, dev_context)
    with open(config_path, "w") as f:
        f.write(job_config)
        info(SUCCESS_MESSAGE % config_path)
    launch_if_open_flagged(config_path, **kwds)
