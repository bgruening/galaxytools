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

from .common_parser import (
    get_common_args,
    HideUnderscoresHelpFormatter,
)

DEFAULT_SLEEP_WAIT = 1
MESSAGE_KEY_NOT_YET_VALID = "[%02d] Provided key not (yet) valid... %s\n"
MESSAGE_INVALID_JSON = "[%02d] No valid json returned... %s\n"
MESSAGE_FETCHING_USER = "[%02d] Connection error fetching user details, exiting with error code. %s\n"
MESSAGE_KEY_NOT_YET_ADMIN = "[%02d] Provided key not (yet) admin... %s\n"
MESSAGE_GALAXY_NOT_YET_UP = "[%02d] Galaxy not up yet... %s\n"
MESSAGE_TIMEOUT = "Failed to contact Galaxy within timeout (%s), exiting with error code.\n"


def _parser():
    """Constructs the parser object"""
    parent = get_common_args(login_required=False)
    parser = ArgumentParser(
        parents=[parent],
        usage="usage: %(prog)s <options>",
        formatter_class=HideUnderscoresHelpFormatter,
        description="Script to sleep and wait for Galaxy to be alive.",
    )
    parser.add_argument(
        "--timeout",
        default=0,
        type=int,
        help="Galaxy startup timeout in seconds. The default value of 0 waits forever",
    )
    parser.add_argument(
        "-a",
        "--api-key",
        "--api_key",
        dest="api_key",
        help="Sleep until key becomes available.",
    )
    parser.add_argument("--ensure-admin", "--ensure_admin", default=False, action="store_true")
    return parser


class SleepCondition:
    def __init__(self):
        self.sleep = True

    def cancel(self):
        self.sleep = False


def galaxy_wait(
    galaxy_url,
    verbose=False,
    timeout=0,
    sleep_condition=None,
    api_key=None,
    ensure_admin=False,
):
    """Pass user_key to ensure it works before returning."""
    if verbose:
        sys.stdout.write(f"calling galaxy_wait with timeout={timeout} ensure_admin={ensure_admin}\n\n\n")
        sys.stdout.flush()

    version_url = galaxy_url + "/api/version"
    if api_key:
        # adding the key to the URL will ensure Galaxy returns invalid responses until
        # the key is available.
        version_url = f"{version_url}?key={api_key}"
        current_user_url = f"{galaxy_url}/api/users/current?key={api_key}"
    else:
        assert not ensure_admin

    if sleep_condition is None:
        sleep_condition = SleepCondition()

    count = 0
    version_obtained = False

    while sleep_condition.sleep:
        try:
            if not version_obtained:
                result = requests.get(version_url)
                if result.status_code == 403:
                    if verbose:
                        sys.stdout.write(MESSAGE_KEY_NOT_YET_VALID % (count, result.__str__()))
                        sys.stdout.flush()
                else:
                    try:
                        result = result.json()
                        if verbose:
                            sys.stdout.write(f"Galaxy Version: {result['version_major']}\n")
                            sys.stdout.flush()
                        version_obtained = True
                    except ValueError:
                        if verbose:
                            sys.stdout.write(MESSAGE_INVALID_JSON % (count, result.__str__()))
                            sys.stdout.flush()

            if version_obtained:
                if ensure_admin:
                    result = requests.get(current_user_url)
                    if result.status_code != 200:
                        if verbose:
                            sys.stdout.write(MESSAGE_FETCHING_USER % (count, result.__str__()))
                            sys.stdout.flush()
                            return False

                    result = result.json()
                    is_admin = result["is_admin"]
                    if is_admin:
                        if verbose:
                            sys.stdout.write("Verified supplied key an admin key.\n")
                            sys.stdout.flush()
                        break
                    else:
                        if verbose:
                            sys.stdout.write(MESSAGE_KEY_NOT_YET_ADMIN % (count, result.__str__()))
                            sys.stdout.flush()
                else:
                    break
        except requests.exceptions.ConnectionError as e:
            if verbose:
                sys.stdout.write(MESSAGE_GALAXY_NOT_YET_UP % (count, unicodify(e)[:100]))
                sys.stdout.flush()
        count += 1

        # If we cannot talk to galaxy and are over the timeout
        if timeout != 0 and count > timeout:
            sys.stderr.write(MESSAGE_TIMEOUT % timeout)
            return False

        time.sleep(DEFAULT_SLEEP_WAIT)

    return True


def main(argv=None):
    """
    Main function
    """
    options = _parser().parse_args(argv)

    galaxy_alive = galaxy_wait(
        galaxy_url=options.galaxy,
        verbose=options.verbose,
        timeout=options.timeout,
        api_key=options.api_key,
        ensure_admin=options.ensure_admin,
    )
    exit_code = 0 if galaxy_alive else 1
    sys.exit(exit_code)


if __name__ == "__main__":
    main()
