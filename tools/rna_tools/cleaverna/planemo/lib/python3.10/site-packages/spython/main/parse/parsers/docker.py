# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json
import os
import re

from spython.logger import bot

from .base import ParserBase


class DockerParser(ParserBase):
    name = "docker"

    def __init__(self, filename="Dockerfile", load=True):
        """a docker parser will read in a Dockerfile and put it into a Recipe
        object.

        Parameters
        ==========
        filename: the Dockerfile to parse. If not defined, defaults to
                  Dockerfile assumed to be in the $PWD.
        load: whether to load the recipe file (default True)

        """
        super(DockerParser, self).__init__(filename, load)

    def parse(self):
        """parse is the base function for parsing the Dockerfile, and extracting
        elements into the correct data structures. Everything is parsed into
        lists or dictionaries that can be assembled again on demand.

        Environment: Since Docker also exports environment as we go,
                     we add environment to the environment section and
                     install

        Labels: include anything that is a LABEL, ARG, or (deprecated)
                maintainer.

        Add/Copy: are treated the same

        """
        parser = None
        previous = None

        for line in self.lines:
            parser = self._get_mapping(line, parser, previous)

            # Parse it, if appropriate
            if parser:
                parser(line)

            previous = line

        # Instantiated by ParserBase
        return self.recipe

    # Setup for each Parser

    def _setup(self, action, line):
        """replace the command name from the group, alert the user of content,
        and clean up empty spaces
        """
        bot.debug("[in]  %s" % line)

        # Replace ACTION at beginning
        line = re.sub("^%s" % action, "", line)

        # Handle continuation lines without ACTION by padding with leading space
        line = " " + line

        # Split into components
        return [x for x in self._split_line(line) if x not in ["", None]]

    # From Parser

    def _from(self, line):
        """get the FROM container image name from a FROM line. If we have
        already seen a FROM statement, this is indicative of adding
        another image (multistage build).

        Parameters
        ==========
        line: the line from the recipe file to parse for FROM
        recipe: the recipe object to populate.
        """
        fromHeader = self._setup("FROM", line)

        # Do we have a multistge build to update the active layer?
        self._multistage(fromHeader[0])

        # Now extract the from header, make args replacements
        self.recipe[self.active_layer].fromHeader = self._replace_from_dict(
            re.sub("AS .+", "", fromHeader[0], flags=re.I), self.args
        )

        if "scratch" in self.recipe[self.active_layer].fromHeader:
            bot.warning("scratch is no longer available on Docker Hub.")
        bot.debug("FROM %s" % self.recipe[self.active_layer].fromHeader)

    # Run and Test Parser

    def _run(self, line):
        """everything from RUN goes into the install list

        Parameters
        ==========
        line: the line from the recipe file to parse for FROM

        """
        line = self._setup("RUN", line)
        self.recipe[self.active_layer].install += line

    def _test(self, line):
        """A healthcheck is generally a test command

        Parameters
        ==========
        line: the line from the recipe file to parse for FROM

        """
        self.recipe[self.active_layer].test = self._setup("HEALTHCHECK", line)

    # Arg Parser

    def _arg(self, line):
        """singularity doesn't have support for ARG, so instead will issue
        a warning to the console for the user to export the variable
        with SINGULARITY prefixed at build.

        Parameters
        ==========
        line: the line from the recipe file to parse for ARG

        """
        line = self._setup("ARG", line)

        # Args are treated like envars, so we add them to install
        environ = self.parse_env([x for x in line if "=" in x])
        self.recipe[self.active_layer].install += environ

        # Try to extract arguments from the line
        for arg in line:
            # An undefined arg cannot be used
            if "=" not in arg:
                bot.warning(
                    "ARG is not supported for Singularity, and must be defined with "
                    "a default to be parsed. Skipping %s" % arg
                )
                continue

            arg, value = arg.split("=", 1)
            arg = arg.strip()
            value = value.strip()
            bot.debug("Updating ARG %s to %s" % (arg, value))
            self.args[arg] = value

    # Env Parser

    def _env(self, line):
        """env will parse a line that beings with ENV, indicative of one or
        more environment variables.

        Parameters
        ==========
        line: the line from the recipe file to parse for ADD

        """
        line = self._setup("ENV", line)

        # Extract environment (list) from the line
        environ = self.parse_env(line)

        # Add to global environment, run during install
        self.recipe[self.active_layer].install += environ

        # Also define for global environment
        self.recipe[self.active_layer].environ += environ

    def parse_env(self, envlist):
        """parse_env will parse a single line (with prefix like ENV removed) to
        a list of commands in the format KEY=VALUE For example:

        ENV PYTHONBUFFER 1 --> [PYTHONBUFFER=1]
        Docker: https://docs.docker.com/engine/reference/builder/#env
        """
        if not isinstance(envlist, list):
            envlist = [envlist]

        exports = []

        for env in envlist:
            pieces = re.split("( |\\\".*?\\\"|'.*?')", env)
            pieces = [p for p in pieces if p.strip()]

            while pieces:
                current = pieces.pop(0)

                if current.endswith("="):
                    # Case 1: ['A='] --> A=
                    nextone = ""

                    # Case 2: ['A=', '"1 2"'] --> A=1 2
                    if pieces:
                        nextone = pieces.pop(0)
                    exports.append("%s%s" % (current, nextone))

                # Case 3: ['A=B']     --> A=B
                elif "=" in current:
                    exports.append(current)

                # Case 4: ENV \\
                elif current.endswith("\\"):
                    continue

                # Case 5: ['A', 'B']  --> A=B
                else:
                    nextone = pieces.pop(0)
                    exports.append("%s=%s" % (current, nextone))

        return exports

    # Add and Copy Parser

    def _copy(self, lines):
        """parse_add will copy multiple files from one location to another.
        This likely will need tweaking, as the files might need to be
        mounted from some location before adding to the image.
        The add command is done for an entire directory. It is also
        possible to have more than one file copied to a destination:
        https://docs.docker.com/engine/reference/builder/#copy
        e.g.: <src> <src> <dest>/
        """
        lines = self._setup("COPY", lines)

        for line in lines:
            # Take into account multistage builds
            layer = None
            if line.startswith("--from"):
                layer = line.strip("--from").split(" ")[0].lstrip("=")
                if layer not in self.recipe:
                    bot.warning(
                        "COPY requested from layer %s, but layer not previously defined."
                        % layer
                    )
                    continue

                # Remove the --from from the line
                line = " ".join([word for word in line.split(" ")[1:] if word])

            values = line.split(" ")
            topath = values.pop()
            for frompath in values:
                self._add_files(frompath, topath, layer)

    def _add(self, lines):
        """Add can also handle https, and compressed files.

        Parameters
        ==========
        line: the line from the recipe file to parse for ADD

        """
        lines = self._setup("ADD", lines)

        for line in lines:
            values = line.split(" ")
            frompath = values.pop(0)

            # Custom parsing for frompath

            # If it's a web address, add to install routine to get
            if frompath.startswith("http"):
                for topath in values:
                    self._parse_http(frompath, topath)

            # Add the file, and decompress in install
            elif re.search("[.](gz|gzip|bz2|xz)$", frompath.strip()):
                for topath in values:
                    self._parse_archive(frompath, topath)

            # Just add the files
            else:
                for topath in values:
                    self._add_files(frompath, topath)

    # File Handling

    def _add_files(self, source, dest, layer=None):
        """add files is the underlying function called to add files to the
        list, whether originally called from the functions to parse archives,
        or https. We make sure that any local references are changed to
        actual file locations before adding to the files list.

        Parameters
        ==========
        source: the source
        dest: the destiation
        """

        # Warn the user Singularity doesn't support expansion
        if "*" in source:
            bot.warning("Singularity doesn't support expansion, * found in %s" % source)

        # Warning if file/folder (src) doesn't exist
        if not os.path.exists(source) and layer is None:
            bot.warning("%s doesn't exist, ensure exists for build" % source)

        # The pair is added to the files as a list
        if not layer:
            self.recipe[self.active_layer].files.append([source, dest])

        # Unless the file is to be copied from a particular layer
        else:
            if layer not in self.recipe[self.active_layer].layer_files:
                self.recipe[self.active_layer].layer_files[layer] = []
            self.recipe[self.active_layer].layer_files[layer].append([source, dest])

    def _parse_http(self, url, dest):
        """will get the filename of an http address, and return a statement
        to download it to some location

        Parameters
        ==========
        url: the source url to retrieve with curl
        dest: the destination folder to put it in the image

        """
        file_name = os.path.basename(url)
        download_path = "%s/%s" % (dest, file_name)
        command = "curl %s -o %s" % (url, download_path)
        self.recipe[self.active_layer].install.append(command)

    def _parse_archive(self, targz, dest):
        """parse_targz will add a line to the install script to extract a
        targz to a location, and also add it to the files.

        Parameters
        ==========
        targz: the targz to extract
        dest: the location to extract it to

        """

        # Add command to extract it
        self.recipe[self.active_layer].install.append("tar -zvf %s %s" % (targz, dest))

        # Ensure added to container files
        return self._add_files(targz, dest)

    # Comments and Default

    def _comment(self, line):
        """Simply add the line to the install as a comment. This function is
        equivalent to default, but added in the case we need future custom
        parsing (meaning a comment is different from a line.

        Parameters
        ==========
        line: the line from the recipe file to parse to INSTALL

        """
        self.recipe[self.active_layer].install.append(line)

    def _default(self, line):
        """the default action assumes a line that is either a command (a
        continuation of a previous, for example) or a comment.

        Parameters
        ==========
        line: the line from the recipe file to parse to INSTALL
        """
        if line.strip().startswith("#"):
            return self._comment(line)
        self.recipe[self.active_layer].install.append(line)

    # Ports and Volumes

    def _volume(self, line):
        """We don't have logic for volume for Singularity, so we add as
        a comment in the install, and a metadata value for the recipe
        object

        Parameters
        ==========
        line: the line from the recipe file to parse to INSTALL

        """
        volumes = self._setup("VOLUME", line)
        if volumes:
            self.recipe[self.active_layer].volumes += volumes
        return self._comment("# %s" % line)

    def _expose(self, line):
        """Again, just add to metadata, and comment in install.

        Parameters
        ==========
        line: the line from the recipe file to parse to INSTALL

        """
        ports = self._setup("EXPOSE", line)
        if ports:
            self.recipe[self.active_layer].ports += ports
        return self._comment("# %s" % line)

    def _stopsignal(self, line):
        """Again, just add to metadata, and comment in install.

        Parameters
        ==========
        line: the line from the recipe file to parse STOPSIGNAL
        """
        return self._comment("# %s" % line)

    # Working Directory

    def _workdir(self, line):
        """A Docker WORKDIR command simply implies to cd to that location

        Parameters
        ==========
        line: the line from the recipe file to parse for WORKDIR

        """
        # Save the last working directory to add to the runscript
        workdir = self._setup("WORKDIR", line)
        workdir_mkdir = "mkdir -p %s" % ("".join(workdir))
        self.recipe[self.active_layer].install.append(workdir_mkdir)
        workdir_cd = "cd %s" % ("".join(workdir))
        self.recipe[self.active_layer].install.append(workdir_cd)
        self.recipe[self.active_layer].workdir = workdir[0]

    # Entrypoint and Command

    def _cmd(self, line):
        """_cmd will parse a Dockerfile CMD command

        eg: CMD /code/run_uwsgi.sh --> /code/run_uwsgi.sh.
            If a list is provided, it's parsed to a list.

        Parameters
        ==========
        line: the line from the recipe file to parse for CMD

        """
        cmd = self._setup("CMD", line)[0]
        self.recipe[self.active_layer].cmd = self._load_list(cmd)

    def _load_list(self, line):
        """load an entrypoint or command, meaning it can be wrapped in a list
        or a regular string. We try loading as json to return an actual
        list. E.g., after _setup, we might go from 'ENTRYPOINT ["one", "two"]'
        to '["one", "two"]', and this function loads as json and returns
        ["one", "two"]
        """
        try:
            line = json.loads(line)
        except Exception:
            pass
        return line

    def _entry(self, line):
        """_entrypoint will parse a Dockerfile ENTRYPOINT command

        Parameters
        ==========
        line: the line from the recipe file to parse for CMD

        """
        entrypoint = self._setup("ENTRYPOINT", line)[0]
        self.recipe[self.active_layer].entrypoint = self._load_list(entrypoint)

    # Labels

    def _label(self, line):
        """_label will parse a Dockerfile label

        Parameters
        ==========
        line: the line from the recipe file to parse for CMD

        """
        label = self._setup("LABEL", line)
        self.recipe[self.active_layer].labels += [label]

    # Main Parsing Functions

    def _get_mapping(self, line, parser=None, previous=None):
        """mapping will take the command from a Dockerfile and return a map
        function to add it to the appropriate place. Any lines that don't
        cleanly map are assumed to be comments.

        Parameters
        ==========
        line: the list that has been parsed into parts with _split_line
        parser: the previously used parser, for context

        Returns
        =======
        function: to map a line to its command group

        """

        # Split the command into cleanly the command and rest
        if not isinstance(line, list):
            line = self._split_line(line)

        # No line we will give function to handle empty line
        if not line:
            return None

        cmd = line[0].upper()

        mapping = {
            "ADD": self._add,
            "ARG": self._arg,
            "COPY": self._copy,
            "CMD": self._cmd,
            "ENTRYPOINT": self._entry,
            "ENV": self._env,
            "EXPOSE": self._expose,
            "FROM": self._from,
            "HEALTHCHECK": self._test,
            "RUN": self._run,
            "WORKDIR": self._workdir,
            "MAINTAINER": self._label,
            "VOLUME": self._volume,
            "LABEL": self._label,
            "STOPSIGNAL": self._stopsignal,
        }

        # If it's a command line, return correct function
        if cmd in mapping:
            return mapping[cmd]

        # If it's a continued line, return previous
        cleaned = self._clean_line(line[-1])
        previous = self._clean_line(previous)

        # if we are continuing from last
        if cleaned.endswith("\\") and parser or previous.endswith("\\"):
            return parser

        return self._default

    def _clean_line(self, line):
        """clean line will remove comments, and strip the line of newlines
        or spaces.

        Parameters
        ==========
        line: the string to parse into parts

        Returns
        =======
        line: a cleaned line

        """
        # A line that is None should return empty string
        line = line or ""
        return line.split("#")[0].strip()
