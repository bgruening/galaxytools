"""test different sets of command line arguments"""

import sys  # for capturing print() output
from collections.abc import Iterator
from contextlib import contextmanager
from io import StringIO

import schema_salad.main as cli_parser


@contextmanager
def captured_output() -> Iterator[tuple[StringIO, StringIO]]:
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


def test_version() -> None:
    args: list[list[str]] = [["--version"], ["-v"]]
    for arg in args:
        with captured_output() as (out, err):
            cli_parser.main(arg)

        response = out.getvalue().strip()  # capture output and strip newline
        assert "Current version" in response


def test_empty_input() -> None:
    # running schema_salad tool without any args
    args: list[str] = []
    with captured_output() as (out, err):
        cli_parser.main(args)

    response = out.getvalue().strip()
    assert "error: too few arguments" in response
