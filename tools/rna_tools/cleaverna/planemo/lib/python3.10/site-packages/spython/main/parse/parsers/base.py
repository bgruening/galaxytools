# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import abc
import os
import re
from copy import deepcopy

from spython.logger import bot
from spython.utils import read_file

from ..recipe import Recipe


class ParserBase:
    """a parser Base is intended to provide helper functions for a parser,
    namely to read lines in files, and otherwise interact with outputs.
    Input should be some recipe (text file to describe a container build)
    and output of parse() is a Recipe (spython.main.parse.recipe.Recipe)
    object, which can be used to write to file, etc.
    """

    def __init__(self, filename, load=True):
        """a generic recipe parser holds the original file, and provides
        shared functions for interacting with files. If the subclass has
        a parse function defined, we parse the filename

        Parameters
        ==========
        filename: the recipe file to parse.
        load: if True, load the filename into the Recipe. If not loaded,
              the user can call self.parse() at a later time.

        """
        self.filename = filename
        self._run_checks()
        self.lines = []

        # Arguments can be used internally, active layer name and number
        self.args = {}
        self.active_layer = "spython-base"
        self.active_layer_num = 1

        # Support multistage builds
        self.recipe = {"spython-base": Recipe(self.filename)}

        if self.filename:
            # Read in the raw lines of the file
            self.lines = read_file(self.filename)

            # If parsing function defined, parse the recipe
            if load:
                self.parse()

    @abc.abstractmethod
    def parse(self):
        """parse is the base function for parsing an input filename, and
        extracting elements into the correct Recipe sections. The exact
        logic and supporting functions will vary based on the recipe type.
        """
        return

    def _run_checks(self):
        """basic sanity checks for the file name (and others if needed) before
        attempting parsing.
        """
        if self.filename is not None:
            # Does the recipe provided exist?
            if not os.path.exists(self.filename):
                bot.exit("Cannot find %s, is the path correct?" % self.filename)

            # Ensure we carry fullpath
            self.filename = os.path.abspath(self.filename)

    # Printing

    def __str__(self):
        """show the user the recipe object, along with the type. E.g.,

        [spython-parser][docker]
        [spython-parser][singularity]

        """
        base = "[spython-parser]"
        if hasattr(self, "name"):
            base = "%s[%s]" % (base, self.name)
        return base

    def __repr__(self):
        return self.__str__()

    # Lines

    def _split_line(self, line):
        """clean a line to prepare it for parsing, meaning separation
        of commands. We remove newlines (from ends) along with extra spaces.

        Parameters
        ==========
        line: the string to parse into parts

        Returns
        =======
        parts: a list of line pieces, the command is likely first

        """
        return [x.strip() for x in line.split(" ", 1)]

    # Multistage

    def _multistage(self, fromHeader):
        """Given a from header, determine if we have a multistage build, and
        update the recipe parser active in case that we do. If we are dealing
        with the first layer and it's named, we also update the default
        name "spython-base" to be what the recipe intended.

        Parameters
        ==========
        fromHeader: the fromHeader parsed from self.from, possibly with AS
        """
        # Derive if there is a named layer
        match = re.search("AS (?P<layer>.+)", fromHeader, flags=re.I)
        if match:
            layer = match.groups("layer")[0].strip()

            # If it's the first layer named incorrectly, we need to rename
            if len(self.recipe) == 1 and list(self.recipe)[0] == "spython-base":
                self.recipe[layer] = deepcopy(self.recipe[self.active_layer])
                del self.recipe[self.active_layer]
            else:
                self.active_layer_num += 1
                self.recipe[layer] = Recipe(self.filename, self.active_layer_num)
            self.active_layer = layer
            bot.debug(
                "Active layer #%s updated to %s"
                % (self.active_layer_num, self.active_layer)
            )

    def _replace_from_dict(self, string, args):
        """Given a lookup of arguments, args, replace any that are found in
        the given string. This is intended to be used to substitute ARGs
        provided in a Dockerfile into other sections, e.g., FROM $BASE

        Parameters
        ==========
        string: an input string to look for replacements
        args: a dictionary to make lookups from

        Returns
        =======
        string: the string with replacements made
        """
        for key, value in args.items():
            if re.search("([$]" + key + "|[$][{]" + key + "[}])", string):
                string = re.sub("([$]" + key + "|[$]{" + key + "[}])", value, string)
        return string
