"""Utilities for scripts in gxformat2."""
from gxformat2.export import from_galaxy_native
from gxformat2.yaml import ordered_load_path


def ensure_format2_from_path(path: str):
    """Load a file from the specified path and ensure it is in format2."""
    return ensure_format2(ordered_load_path(path))


def ensure_format2(workflow_dict: dict, ensure_labels: bool = False):
    """Consume a dictionary and ensure the result is format2.

    So convert from ga if needed.
    """
    if workflow_dict.get("a_galaxy_workflow") == "true":
        workflow_dict = from_galaxy_native(workflow_dict)
    return workflow_dict
