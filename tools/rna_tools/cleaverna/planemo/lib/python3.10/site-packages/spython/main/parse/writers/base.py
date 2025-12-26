# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import os
import tempfile

from spython.logger import bot
from spython.utils import write_file


class WriterBase:
    def __init__(self, recipe=None):
        """a writer base will take a recipe object (parser.base.Recipe) and
        provide helpers for writing to file.

        Parameters
        ==========
        recipe: the recipe instance to parse

        """
        self.recipe = recipe

    def write(self, output_file=None, force=False):
        """convert a recipe to a specified format, and write to file, meaning
        we use the loaded recipe to write to an output file.
        If the output file is not specified, a temporary file is used.

        Parameters
        ==========
        output_file: the file to save to, not required (estimates default)
        force: if True, if file exists, over-write existing file

        """
        if output_file is None:
            output_file = self._get_conversion_outfile()

        # Cut out early if file exists and we aren't overwriting
        if os.path.exists(output_file) and not force:
            bot.exit("%s exists, and force is False." % output_file)

        # Do the conversion if function is provided by subclass
        if hasattr(self, "convert"):
            converted = self.convert()
            bot.info("Saving to %s" % output_file)
            write_file(output_file, converted)

    def _get_conversion_outfile(self):
        """a helper function to return a conversion temporary output file
        based on kind of conversion

        Parameters
        ==========
        convert_to: a string either docker or singularity, if a different

        """
        prefix = "spythonRecipe"
        if hasattr(self, "name"):
            prefix = self.name
        suffix = next(tempfile._get_candidate_names())
        return "%s.%s" % (prefix, suffix)

    # Printing

    def __str__(self):
        """show the user the recipe object, along with the type. E.g.,

        [spython-writer][docker]
        [spython-writer][singularity]

        """
        base = "[spython-writer]"
        if hasattr(self, "name"):
            base = "%s[%s]" % (base, self.name)
        return base

    def __repr__(self):
        return self.__str__()
