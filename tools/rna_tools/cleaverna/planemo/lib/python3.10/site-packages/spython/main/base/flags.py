# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


"""
GLOBAL OPTIONS:
    -d|--debug    Print debugging information
    -h|--help     Display usage summary
    -s|--silent   Only print errors
    -q|--quiet    Suppress all normal output
       --version  Show application version
    -v|--verbose  Increase verbosity +1
    -x|--sh-debug Print shell wrapper debugging information

GENERAL COMMANDS:
    help       Show additional help for a command or container
    selftest   Run some self tests for singularity install

CONTAINER USAGE COMMANDS:
    exec       Execute a command within container
    run        Launch a runscript within container
    shell      Run a Bourne shell within container
    test       Launch a testscript within container

CONTAINER MANAGEMENT COMMANDS:
    apps       List available apps within a container
    bootstrap  *Deprecated* use build instead
    build      Build a new Singularity container
    check      Perform container lint checks
    inspect    Display container's metadata
    mount      Mount a Singularity container image
    pull       Pull a Singularity/Docker container to $PWD
    siflist    list data object descriptors of a SIF container image
    sign       Sign a group of data objects in container
    verify     Verify the crypto signature of group of data objects in container

COMMAND GROUPS:
    capability User's capabilities management command group
    image      Container image command group
    instance   Persistent instance command group

"""


def parse_verbosity(self, args):
    """parse_verbosity will take an argument object, and return the args
    passed (from a dictionary) to a list

    Parameters
    ==========
    args: the argparse argument objects

    """

    flags = []

    if args.silent:
        flags.append("--silent")
    elif args.quiet:
        flags.append("--quiet")
    elif args.debug:
        flags.append("--debug")
    elif args.verbose:
        flags.append("-" + "v" * args.verbose)

    return flags
