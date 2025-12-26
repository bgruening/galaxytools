# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def helpcmd(self, command=None):
    """help prints the general function help, or help for a specific command

    Parameters
    ==========
    command: the command to get help for, if none, prints general help

    """
    from spython.utils import check_install

    check_install()

    cmd = ["singularity", "--help"]
    if command is not None:
        cmd.append(command)
    return self._run_command(cmd)
