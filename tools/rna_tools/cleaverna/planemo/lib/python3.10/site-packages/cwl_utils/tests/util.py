import atexit
import os
import shutil
from contextlib import ExitStack
from importlib.resources import as_file, files
from pathlib import Path

import pytest


def get_path(filename: str) -> Path:
    """Get the filepath for a given test file."""
    # normalizing path depending on OS or else it will cause problem when joining path
    filename = os.path.normpath(filename)
    filepath = None
    try:
        file_manager = ExitStack()
        atexit.register(file_manager.close)
        traversable = files("cwl-utils") / filename
        filepath = file_manager.enter_context(as_file(traversable))
    except ModuleNotFoundError:
        pass
    if not filepath or not filepath.is_file():
        filepath = Path(os.path.dirname(__file__), os.pardir, filename)
    return filepath.resolve()


def get_data(filename: str) -> str:
    return str(get_path(filename))


needs_docker = pytest.mark.skipif(
    not bool(shutil.which("docker")),
    reason="Requires the docker executable on the system path.",
)

needs_singularity = pytest.mark.skipif(
    not bool(shutil.which("singularity")),
    reason="Requires the singularity executable on the system path.",
)

needs_podman = pytest.mark.skipif(
    not bool(shutil.which("podman")),
    reason="Requires the podman executable on the system path.",
)

needs_udocker = pytest.mark.skipif(
    not bool(shutil.which("udocker")),
    reason="Requires the udocker executable on the system path.",
)
