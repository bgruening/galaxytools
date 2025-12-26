#!/usr/bin/env python
"""Utility to do a blocking sleep until a Galaxy instance is responsive.
This is useful in docker images, in RUN steps, where one needs to wait
for a currently starting Galaxy to be alive, before API requests can be
made successfully.
The script functions by making repeated requests to
``http(s)://fqdn/api/version``, an API which requires no authentication
to access."""

import sys
import time
from argparse import ArgumentParser

import requests
from galaxy.util import unicodify

try:
    from .common_parser import get_common_args
except ImportError:
    # This won't stay in Planemo long, main no longer functional.
    get_common_args = None

DEFAULT_SLEEP_WAIT = 1


def _parser():
    """Constructs the parser object"""
    parent = get_common_args(login_required=False)
    parser = ArgumentParser(
        parents=[parent],
        usage="usage: python %(prog)s <options>",
        description="Script to sleep and wait for Galaxy to be alive.",
    )
    parser.add_argument(
        "--timeout", default=0, type=int, help="Galaxy startup timeout in seconds. The default value of 0 waits forever"
    )
    return parser


def _parse_cli_options():
    """
    Parse command line options, returning `parse_args` from `ArgumentParser`.
    """
    parser = _parser()
    return parser.parse_args()


class SleepCondition:
    def __init__(self):
        self.sleep = True

    def cancel(self):
        self.sleep = False


def sleep(galaxy_url, verbose=False, timeout=0, sleep_condition=None):
    if sleep_condition is None:
        sleep_condition = SleepCondition()

    count = 0
    while sleep_condition.sleep:
        try:
            result = requests.get(galaxy_url + "/api/version")
            try:
                result = result.json()
                if verbose:
                    sys.stdout.write("Galaxy Version: %s\n" % result["version_major"])
                    sys.stdout.flush()
                break
            except ValueError:
                if verbose:
                    sys.stdout.write("[%02d] No valid json returned... %s\n" % (count, result.__str__()))
                    sys.stdout.flush()
        except requests.exceptions.ConnectionError as e:
            if verbose:
                sys.stdout.write("[%02d] Galaxy not up yet... %s\n" % (count, unicodify(e)[:100]))
                sys.stdout.flush()
        count += 1

        # If we cannot talk to galaxy and are over the timeout
        if timeout != 0 and count > timeout:
            sys.stderr.write("Failed to contact Galaxy\n")
            return False

        time.sleep(DEFAULT_SLEEP_WAIT)

    return True


def main():
    """
    Main function
    """
    options = _parse_cli_options()

    galaxy_alive = sleep(galaxy_url=options.galaxy, verbose=options.verbose, timeout=options.timeout)
    exit_code = 0 if galaxy_alive else 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
