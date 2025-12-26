"""Utilities for using virtualenv as library and planemo command."""

import os
import sys
from typing import Optional

from galaxy.util.commands import which

DEFAULT_PYTHON_VERSION = os.environ.get("PLANEMO_DEFAULT_PYTHON_VERSION", "3")


def create_command(virtualenv_path: str, galaxy_python_version: Optional[str] = None) -> str:
    """If virtualenv is on Planemo's path use it, otherwise use the planemo
    subcommand virtualenv to create the virtualenv.
    """
    # Create a virtualenv with the selected python version.
    if galaxy_python_version is None:
        galaxy_python_version = DEFAULT_PYTHON_VERSION
    python = which("python%s" % galaxy_python_version)
    if python:
        python = os.path.abspath(python)
    else:
        python = sys.executable or "python"
    virtualenv_on_path = which("virtualenv")
    if virtualenv_on_path:
        command = [virtualenv_on_path, virtualenv_path, "-p", python]
    else:
        command = [python, "-m", "venv", virtualenv_path]
    return " ".join(command)
