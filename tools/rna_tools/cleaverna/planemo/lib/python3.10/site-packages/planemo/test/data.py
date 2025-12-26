"""Utilities related to reasoning about test data."""

import os


def find_test_data_directory(tool_paths, **kwds):
    path = "."
    if len(tool_paths) > 0:
        path = tool_paths[0]

    # Find test data directory associated with path.
    test_data = kwds.get("test_data", None)
    if test_data:
        return os.path.abspath(test_data)
    else:
        test_data = search_tool_path_for(path, "test-data")
        if test_data:
            return test_data


def search_tool_path_for(path, target, extra_paths=[]):
    """Check for presence of a target in different artifact directories."""
    if not os.path.isdir(path):
        tool_dir = os.path.dirname(path)
    else:
        tool_dir = path
    possible_dirs = [tool_dir, "."] + extra_paths
    for possible_dir in possible_dirs:
        possible_path = os.path.join(possible_dir, target)
        if os.path.exists(possible_path):
            return os.path.abspath(possible_path)
    return None


__all__ = ("search_tool_path_for",)
