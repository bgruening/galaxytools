# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os

from spython.image import ImageBase
from spython.utils.fileio import read_file


class Instance(ImageBase):
    def __init__(self, image, start=True, name=None, quiet=True, **kwargs):
        """An instance is an image running as an instance with services.
        This class has functions appended under cmd/__init__ and is
        instantiated when the user calls Client.

        Parameters
        ==========
        image: the Singularity image uri to parse (required)
        start: boolean to start the instance (default is True)
        name: a name for the instance (will generate RobotName
                if not provided)
        """
        super(Instance, self).__init__()
        self.parse_image_name(image)
        self.generate_name(name)

        # Update metadats from arguments
        self._update_metadata(kwargs)
        self.options = kwargs.get("options", [])
        self.args = kwargs.get("args", [])
        self.cmd = []
        self.quiet = quiet

        # Start the instance
        if start:
            self.start(quiet=quiet, **kwargs)

    # Unique resource identifier

    def generate_name(self, name=None):
        """generate a Robot Name for the instance to use, if the user doesn't
        supply one.
        """
        # If no name provided, use robot name
        if name is None:
            name = self.RobotNamer.generate()
        # dash allowed in instance name.
        # authorizedChars = `^[a-zA-Z0-9._-]+$` from instance_linux.go
        # self.name = name.replace("-", "_")
        self.name = name

    def parse_image_name(self, image):
        """
        simply split the uri from the image. Singularity handles
        parsing of registry, namespace, image.

        Parameters
        ==========
        image: the complete image uri to load (e.g., docker://ubuntu)

        """
        self._image = image
        self.protocol = "instance"

    def get_uri(self):
        """return the image uri (instance://) along with it's name"""
        return self.__str__()

    # Logs

    def read_logs_out(self):
        """Read output log file, if it exists"""
        if hasattr(self, "log_out_path"):
            if os.path.exists(self.log_out_path):
                return read_file(self.log_out_path)

    def read_logs_err(self):
        """Read error log file, if it exists"""
        if hasattr(self, "log_err_path"):
            if os.path.exists(self.log_err_path):
                return read_file(self.log_err_path)

    # Metadata

    def _update_metadata(self, kwargs=None):
        """Extract any additional attributes to hold with the instance
        from kwargs
        """

        # If not given metadata, use instance.list to get it for container
        if kwargs is None and hasattr(self, "name"):
            kwargs = self._list(self.name, quiet=True, return_json=True)

        # Add acceptable arguments
        for arg in ["pid", "name", "ip_address", "log_err_path", "log_out_path", "img"]:
            # Skip over non-iterables:
            if arg in kwargs:
                setattr(self, arg, kwargs[arg])

        if "image" in kwargs:
            self._image = kwargs["image"]
        elif "container_image" in kwargs:
            self._image = kwargs["container_image"]

    def __str__(self):
        if hasattr(self, "name"):
            if self.protocol:
                return "%s://%s" % (self.protocol, self.name)
        return os.path.basename(self._image)

    def __repr__(self):
        return self.__str__()
