#!/usr/bin/env python

import argparse
import os

DEFAULT_JOB_SLEEP = 3


class HideUnderscoresHelpFormatter(argparse.HelpFormatter):
    def add_arguments(self, actions):
        for action in actions:
            action.option_strings = list(s for s in action.option_strings if "_" not in s)
            self.add_argument(action)


class RawDescriptionHideUnderscoresHelpFormatter(HideUnderscoresHelpFormatter, argparse.RawDescriptionHelpFormatter):
    pass


class ArgumentDefaultsHideUnderscoresHelpFormatter(
    HideUnderscoresHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    pass


def add_verbosity_argument(parser_or_group):
    parser_or_group.add_argument("-v", "--verbose", help="Increase output verbosity.", action="store_true")


def add_log_file_argument(parser_or_group):
    parser_or_group.add_argument(
        "--log-file",
        "--log_file",
        dest="log_file",
        help="Where the log file should be stored. " "Default is a file in your system's temp folder",
        default=None,
    )


def get_common_args(login_required=True, log_file=False):
    parser = argparse.ArgumentParser(add_help=False)
    general_group = parser.add_argument_group("General options")
    add_verbosity_argument(general_group)
    if log_file:
        add_log_file_argument(general_group)

    con_group = parser.add_argument_group("Galaxy connection")
    default_galaxy = os.environ.get("EPHEMERIS_GALAXY") or "http://localhost:8080"
    con_group.add_argument(
        "-g",
        "--galaxy",
        help="Target Galaxy instance URL/IP address",
        default=default_galaxy,
    )

    if login_required:
        con_group.add_argument("-u", "--user", help="Galaxy user email address")
        con_group.add_argument("-p", "--password", help="Password for the Galaxy user")
        con_group.add_argument(
            "-a",
            "--api-key",
            "--api_key",
            dest="api_key",
            help="Galaxy admin user API key (required if not defined in the tools list file)",
        )

    return parser
