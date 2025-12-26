"""The module describes a CLI framework extending ``click``."""

import functools
import os
import sys
from typing import (
    Any,
    Callable,
    Dict,
    List,
)

import click
from click.core import (
    Command,
    Context,
)

from planemo import __version__
from planemo.context import (
    configure_standard_planemo_logging,
    PlanemoContext,
)
from planemo.exit_codes import ExitCodeException
from planemo.galaxy import profiles
from .config import OptionSource

CONTEXT_SETTINGS = dict(auto_envvar_prefix="PLANEMO")
COMMAND_ALIASES = {
    "l": "lint",
    "o": "open",
    "t": "test",
    "s": "serve",
}


class PlanemoCliContext(PlanemoContext):
    """Describe context of Planemo CLI computation.

    Extend PlanemoContext with operations for CLI concerns (exit code) and
    for interacting with the click library.
    """

    def _log_message(self, message: str) -> None:
        click.echo(message, file=sys.stderr)

    def exit(self, exit_code: int):
        """Exit planemo with the supplied exit code."""
        self.vlog("Exiting planemo with exit code [%d]" % exit_code)
        raise ExitCodeException(exit_code)


pass_context = click.make_pass_decorator(PlanemoCliContext, ensure=True)
cmd_folder = os.path.abspath(os.path.join(os.path.dirname(__file__), "commands"))


def list_cmds() -> List[str]:
    """List planemo commands from commands folder."""
    rv = []
    for filename in os.listdir(cmd_folder):
        if filename.endswith(".py") and filename.startswith("cmd_"):
            rv.append(filename[len("cmd_") : -len(".py")])
    rv.sort()
    return rv


def name_to_command(name: str) -> Command:
    """Convert a subcommand name to the cli function for that command.

    Command <X> is defined by the method 'planemo.commands.cmd_<x>:cli',
    this method uses `__import__` to load and return that method.
    """
    try:
        mod_name = "planemo.commands.cmd_" + name
        mod = __import__(mod_name, None, None, ["cli"])
    except ImportError as e:
        # Check if this is an unknown command vs. a real import error
        available_commands = list_cmds()
        if name not in available_commands:
            # User tried to run a command that doesn't exist
            error_msg = f"No such command '{name}'."

            # Try to find similar commands to suggest
            similar_commands = _find_similar_commands(name, available_commands)
            if similar_commands:
                error_msg += "\n\nDid you mean one of these?\n"
                for cmd in similar_commands[:5]:  # Show at most 5 suggestions
                    error_msg += f"  - {cmd}\n"

            error_msg += "\nRun 'planemo --help' to see all available commands."
            raise click.UsageError(error_msg)
        else:
            # This is a real import error for a command that should exist
            raise Exception(f"Problem loading command {name}, exception {e}")
    return mod.cli


def _find_similar_commands(name: str, available_commands: List[str]) -> List[str]:
    """Find commands similar to the given name using simple string matching."""
    similar = []

    # First, find commands that start with the same prefix
    for cmd in available_commands:
        if cmd.startswith(name[:3]) or name.startswith(cmd[:3]):
            similar.append(cmd)

    # If we didn't find any prefix matches, look for substring matches
    if not similar:
        for cmd in available_commands:
            if name in cmd or cmd in name:
                similar.append(cmd)

    # If still no matches, look for commands with similar words
    if not similar:
        name_parts = set(name.split("_"))
        for cmd in available_commands:
            cmd_parts = set(cmd.split("_"))
            if name_parts & cmd_parts:  # Check if there's any word overlap
                similar.append(cmd)

    return similar


class PlanemoCLI(click.MultiCommand):
    def list_commands(self, ctx: Context) -> List[str]:
        return list_cmds()

    def get_command(self, ctx: Context, name: str) -> Command:
        if name in COMMAND_ALIASES:
            name = COMMAND_ALIASES[name]
        return name_to_command(name)


def command_function(f: Callable) -> Callable:
    """Extension point for processing kwds after click callbacks."""

    @functools.wraps(f)
    def handle_blended_options(*args, **kwds):
        profile = kwds.get("profile", None)
        if profile:
            ctx = args[0]
            profile_defaults = profiles.ensure_profile(ctx, profile, **kwds)
            _setup_profile_options(ctx, profile_defaults, kwds)

        try:
            return f(*args, **kwds)
        except ExitCodeException as e:
            sys.exit(e.exit_code)

    return pass_context(handle_blended_options)


def _setup_profile_options(ctx: PlanemoCliContext, profile_defaults: Dict[str, Any], kwds: Dict[str, Any]) -> None:
    for key, value in profile_defaults.items():
        option_present = key in kwds
        option_cli_specified = option_present and (ctx.get_option_source(key) == OptionSource.cli)
        use_profile_option = not option_present or not option_cli_specified
        if use_profile_option:
            kwds[key] = value
            ctx.set_option_source(key, OptionSource.profile, force=True)


@click.command(cls=PlanemoCLI, context_settings=CONTEXT_SETTINGS)
@click.version_option(__version__)
@click.option("-v", "--verbose", is_flag=True, help="Enables verbose mode.")
@click.option(
    "--config", default="~/.planemo.yml", envvar="PLANEMO_GLOBAL_CONFIG_PATH", help="Planemo configuration YAML file."
)
@click.option("--directory", default="~/.planemo", envvar="PLANEMO_GLOBAL_WORKSPACE", help="Workspace for planemo.")
@pass_context
def planemo(ctx, config, directory, verbose, configure_logging=True):
    """A command-line toolkit for building tools and workflows for Galaxy.

    Check out the full documentation for Planemo online
    http://planemo.readthedocs.org or open with ``planemo docs``.

    All the individual planemo commands support the ``--help`` option, for
    example use ``planemo lint --help`` for more details on checking tools.
    """
    ctx.verbose = verbose
    if configure_logging:
        configure_standard_planemo_logging(verbose)
    ctx.planemo_config = os.path.expanduser(config)
    ctx.planemo_directory = os.path.expanduser(directory)


__all__ = (
    "command_function",
    "list_cmds",
    "name_to_command",
    "planemo",
    "PlanemoCliContext",
)
