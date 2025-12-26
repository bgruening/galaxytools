# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import re

from spython.logger import bot

from .base import WriterBase


class SingularityWriter(WriterBase):
    name = "singularity"

    def __init__(self, recipe=None):  # pylint: disable=useless-super-delegation
        """a SingularityWriter will take a Recipe as input, and write
        to a Singularity recipe file.

        Parameters
        ==========
        recipe: the Recipe object to write to file.

        """
        super(SingularityWriter, self).__init__(recipe)

    def validate(self):
        """validate that all (required) fields are included for the Docker
        recipe. We minimimally just need a FROM image, and must ensure
        it's in a valid format. If anything is missing, we exit with error.
        """
        if self.recipe is None:
            bot.exit("Please provide a Recipe() to the writer first.")

    def convert(self, runscript="/bin/bash", force=False):
        """docker2singularity will return a Singularity build recipe based on
        a the loaded recipe object. It doesn't take any arguments as the
        recipe object contains the sections, and the calling function
        determines saving / output logic.
        """
        self.validate()

        # Write single recipe that includes all layer
        recipe = []

        # Number of layers
        num_layers = len(self.recipe)
        count = 0

        # Write each layer to new file
        for stage, parser in self.recipe.items():
            # Set the first and active stage
            self.stage = stage

            # From header is required
            if parser.fromHeader is None:
                bot.exit("Singularity recipe requires a from header.")

            recipe += ["\n\n\nBootstrap: docker"]
            recipe += ["From: %s" % parser.fromHeader]
            recipe += ["Stage: %s\n\n\n" % stage]

            # TODO: stopped here - bug with files being found
            # Add global files, and then layer files
            recipe += self._create_section("files")
            for layer, files in parser.layer_files.items():
                recipe += create_keyval_section(files, "files", layer)

            # Sections with key value pairs
            recipe += self._create_section("labels")
            recipe += self._create_section("install", "post")
            recipe += self._create_section("environ", "environment")

            # If we are at the last layer, write the runscript
            if count == num_layers - 1:
                runscript = self._create_runscript(runscript, force)

                # If a working directory was used, add it as a cd
                if parser.workdir is not None:
                    runscript = ["cd " + parser.workdir] + [runscript]

                # Finish the recipe, also add as startscript
                recipe += finish_section(runscript, "runscript")
                recipe += finish_section(runscript, "startscript")

                if parser.test is not None:
                    recipe += finish_section(parser.test, "test")
            count += 1

        # Clean up extra white spaces
        recipe = "\n".join(recipe).replace("\n\n", "\n").strip("\n")
        return recipe.rstrip()

    def _create_runscript(self, default="/bin/bash", force=False):
        """create_entrypoint is intended to create a singularity runscript
        based on a Docker entrypoint or command. We first use the Docker
        ENTRYPOINT, if defined. If not, we use the CMD. If neither is found,
        we use function default.

        Parameters
        ==========
        default: set a default entrypoint, if the container does not have
                 an entrypoint or cmd.
        force: If true, use default and ignore Dockerfile settings
        """
        entrypoint = default

        # Only look at Docker if not enforcing default
        if not force:
            if self.recipe[self.stage].entrypoint is not None:
                # The provided entrypoint can be a string or a list
                if isinstance(self.recipe[self.stage].entrypoint, list):
                    entrypoint = " ".join(self.recipe[self.stage].entrypoint)
                else:
                    entrypoint = "".join(self.recipe[self.stage].entrypoint)

            if self.recipe[self.stage].cmd is not None:
                if isinstance(self.recipe[self.stage].cmd, list):
                    entrypoint = (
                        entrypoint + " " + " ".join(self.recipe[self.stage].cmd)
                    )
                else:
                    entrypoint = entrypoint + " " + "".join(self.recipe[self.stage].cmd)

        # Entrypoint should use exec
        if not entrypoint.startswith("exec"):
            entrypoint = "exec %s" % entrypoint

        # Should take input arguments into account
        if not re.search('"?[$]@"?', entrypoint):
            entrypoint = '%s "$@"' % entrypoint
        return entrypoint

    def _create_section(self, attribute, name=None, stage=None):
        """create a section based on key, value recipe pairs,
         This is used for files or label

        Parameters
        ==========
        attribute: the name of the data section, either labels or files
        name: the name to write to the recipe file (e.g., %name).
              if not defined, the attribute name is used.

        """

        # Default section name is the same as attribute
        if name is None:
            name = attribute

        # Put a space between sections
        section = ["\n"]

        # Only continue if we have the section and it's not empty
        try:
            section = getattr(self.recipe[self.stage], attribute)
        except AttributeError:
            bot.debug("Recipe does not have section for %s" % attribute)
            return section

        # if the section is empty, don't print it
        if not section:
            return section

        # Files
        if attribute in ["files", "labels"]:
            return create_keyval_section(section, name, stage)

        # An environment section needs exports
        if attribute in ["environ"]:
            return create_env_section(section, name)

        # Post, Setup
        return finish_section(section, name)


def finish_section(section, name):
    """finish_section will add the header to a section, to finish the recipe
    take a custom command or list and return a section.

    Parameters
    ==========
    section: the section content, without a header
    name: the name of the section for the header

    """

    if not isinstance(section, list):
        section = [section]

    # Convert USER lines to change user
    lines = []
    for line in section:
        if re.search("^USER", line):
            username = line.replace("USER", "", 1).rstrip()
            line = "su - %s" % username + " # " + line
        lines.append(line)

    header = ["%" + name]
    return header + lines


def create_keyval_section(pairs, name, layer):
    """create a section based on key, value recipe pairs,
     This is used for files or label

    Parameters
    ==========
    section: the list of values to return as a parsed list of lines
    name: the name of the section to write (e.g., files)
    layer: if a layer name is provided, name section
    """
    if layer:
        section = ["%" + name + " from %s" % layer]
    else:
        section = ["%" + name]
    for pair in pairs:
        section.append(" ".join(pair).strip().strip("\\"))
    return section


def create_env_section(pairs, name):
    """environment key value pairs need to be joined by an equal, and
     exported at the end.

    Parameters
    ==========
    section: the list of values to return as a parsed list of lines
    name: the name of the section to write (e.g., files)

    """
    section = ["%" + name]
    for pair in pairs:
        section.append("export %s" % pair)
    return section
