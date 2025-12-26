# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


from spython.logger import bot
from spython.utils import check_install, get_singularity_version

from .command import generate_bind_list, init_command, run_command
from .flags import parse_verbosity
from .generate import RobotNamer
from .logger import init_level, println
from .sutils import get_filename, get_uri, load, setenv


class Client:
    def __str__(self):
        base = "[singularity-python]"
        if hasattr(self, "simage"):
            if self.simage.image not in [None, ""]:
                base = "%s[%s]" % (base, self.simage)
        return base

    def __repr__(self):
        return self.__str__()

    def __init__(self):
        """
        The base client for singularity, will have commands added to it.
        upon init, store verbosity requested in environment MESSAGELEVEL.
        """
        self._init_level()

    def version(self):
        """
        Shortcut to get_singularity_version, takes no arguments.
        """
        return get_singularity_version()

    def _check_install(self):
        """
        Ensure that singularity is installed, and exit if not.
        """
        if check_install() is not True:
            bot.exit("Cannot find Singularity! Is it installed?")


# Image Utils
Client.load = load
Client._get_filename = get_filename
Client._get_uri = get_uri
Client.setenv = setenv

# Commands
Client._generate_bind_list = generate_bind_list
Client._init_command = init_command
Client._run_command = run_command

# Flags and Logger
Client._parse_verbosity = parse_verbosity
Client._println = println
Client._init_level = init_level
Client.RobotNamer = RobotNamer()
