#!/usr/bin/env python
"""Wrapper around pytest to execute the bioblend Galaxy test suite against fixed instance.

By default all Galaxy tests will run but a smaller subset can be executed by setting
the environment variable ``BIOBLEND_TEST_SUITE`` to ``quick``.
"""
import os
import sys
from typing import (
    NoReturn,
    Optional,
)

try:
    import pytest
except ImportError:
    pytest = None

DIRECTORY = os.path.abspath(os.path.dirname(__file__))
BIOBLEND_TEST_SUITE = os.environ.get("BIOBLEND_TEST_SUITE", "full")

quick_tests = [
    "TestGalaxyRoles.py",
    "TestGalaxyRoles.py",
    "TestGalaxyUsers.py",
    "TestGalaxyToolData.py",
    "TestGalaxyTools.py::TestGalaxyTools::test_get_tools",  # Test single upload command.
]


def main(args: Optional[list[str]] = None) -> NoReturn:
    """Entry point that delegates to pytest.main."""
    if pytest is None:
        raise Exception("pytest is required to use this script.")
    if args is None:
        args = sys.argv[1:]
    if len(args) < 2:
        if BIOBLEND_TEST_SUITE == "full":
            args.append(os.path.join(DIRECTORY))
        else:
            for quick_test in quick_tests:
                args.append(os.path.join(DIRECTORY, quick_test))
    sys.exit(pytest.main(args))


if __name__ == "__main__":
    main()
