# SPDX-License-Identifier: Apache-2.0
"""Test __meta__ properties."""

from cwl_utils.__meta__ import __version__


def test_graph_split() -> None:
    """Confirm that __version__ exists and is a string."""
    assert __version__
    assert isinstance(__version__, str)
