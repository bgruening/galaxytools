"""Utilities for calling Galaxy scripts."""

import os
import shlex
import string
from typing import (
    Any,
    Dict,
    Optional,
    TYPE_CHECKING,
)

from galaxy.util.commands import shell

from planemo.io import (
    info,
    shell_join,
)
from planemo.virtualenv import (
    create_command,
    DEFAULT_PYTHON_VERSION,
)

if TYPE_CHECKING:
    from planemo.galaxy.config import LocalGalaxyConfig

# Activate galaxy's virtualenv if present (needed for tests say but not for
# server because run.sh does this).
ACTIVATE_COMMAND = (
    'if [ -e "$GALAXY_VIRTUAL_ENV" ]; then . "$GALAXY_VIRTUAL_ENV"/bin/activate; '
    'echo "Activated a virtualenv for Galaxy"; echo "$VIRTUAL_ENV"; '
    'else echo "Failed to activate virtualenv."; fi'
)
CREATE_COMMAND_TEMPLATE = string.Template(
    'if [ ! -e "$GALAXY_VIRTUAL_ENV" ]; then $create_virtualenv; echo "Created virtualenv"; fi',
)
PRINT_VENV_COMMAND = shell_join(
    r'echo "Set \$GALAXY_VIRTUAL_ENV to $GALAXY_VIRTUAL_ENV"',
    (
        'if [ -e "$GALAXY_VIRTUAL_ENV" ]; '
        'then echo "Virtual environment directory exists."; '
        'else echo "Virtual environment directory does not exist."; fi'
    ),
)


CACHED_VIRTUAL_ENV_COMMAND = "if [ -d .venv ]; then GALAXY_VIRTUAL_ENV=.venv; else GALAXY_VIRTUAL_ENV=%s; fi"
UNCACHED_VIRTUAL_ENV_COMMAND = "GALAXY_VIRTUAL_ENV=.venv"


def setup_venv(ctx, kwds: Dict[str, Any], config: Optional["LocalGalaxyConfig"] = None):
    if kwds.get("skip_venv", False):
        return ""

    create_template_params = {
        "create_virtualenv": create_command("$GALAXY_VIRTUAL_ENV", kwds.get("galaxy_python_version"))
    }
    return shell_join(
        locate_galaxy_virtualenv(ctx, kwds, config),
        PRINT_VENV_COMMAND if ctx.verbose else None,
        CREATE_COMMAND_TEMPLATE.safe_substitute(create_template_params),
        PRINT_VENV_COMMAND if ctx.verbose else None,
        ACTIVATE_COMMAND,
        "bash -c 'pwd'" if ctx.verbose else None,
        "bash -c 'which python'" if ctx.verbose else None,
        "bash -c 'which pip'" if ctx.verbose else None,
        "bash -c 'echo $GALAXY_VIRTUAL_ENV'" if ctx.verbose else None,
        "bash -c 'echo $VIRTUAL_ENV'" if ctx.verbose else None,
        "bash -c 'ls -a'" if ctx.verbose else None,
    )


def locate_galaxy_virtualenv(ctx, kwds: Dict[str, Any], config: Optional["LocalGalaxyConfig"] = None):
    virtual_env_locs = []
    if os.environ.get("GALAXY_VIRTUAL_ENV"):
        venv_command = ""
        virtual_env_locs.append(os.environ["GALAXY_VIRTUAL_ENV"])
    elif not kwds.get("no_cache_galaxy", False):
        workspace = ctx.workspace
        galaxy_branch = kwds.get("galaxy_branch") or "master"
        shared_venv_path = os.path.join(workspace, "gx_venv")
        galaxy_python_version = kwds.get("galaxy_python_version") or DEFAULT_PYTHON_VERSION
        shared_venv_path = f"{shared_venv_path}_{galaxy_python_version}"
        if galaxy_branch != "master":
            shared_venv_path = f"{shared_venv_path}_{galaxy_branch}"

        virtual_env_locs += [".venv", shared_venv_path]
        venv_command = CACHED_VIRTUAL_ENV_COMMAND % shlex.quote(shared_venv_path)
    else:
        virtual_env_locs.append(".venv")
        venv_command = UNCACHED_VIRTUAL_ENV_COMMAND
    if config:
        config._virtual_env_locs = virtual_env_locs
    return shell_join(
        venv_command,
        "export GALAXY_VIRTUAL_ENV",
    )


def run_galaxy_command(ctx, command, env, action):
    """Run Galaxy command with informative verbose logging."""
    message = f"{action} with command [{command}]"
    # info not working in pytest+Github actions the way it did in nose?
    info(message)
    ctx.vlog("With environment variables:")
    ctx.vlog("============================")
    for key, value in env.items():
        ctx.vlog(f'{key}="{value}"')
    ctx.vlog("============================")
    exit_code = shell(command, env=env)
    ctx.vlog("run command exited with return code %s" % exit_code)
    return exit_code


__all__ = (
    "setup_venv",
    "run_galaxy_command",
)
