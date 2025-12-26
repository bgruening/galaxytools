# Singularity Image utils for interacting with the Image/Instance
#           classes from the client

# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import os
import re

from spython.logger import bot


def load(self, image=None):
    """load an image, either an actual path on the filesystem or a uri.

    Parameters
    ==========
    image: the image path or uri to load (e.g., docker://ubuntu

    """
    from spython.image import Image
    from spython.instance import Instance

    self.simage = Image(image)

    if image is not None:
        if image.startswith("instance://"):
            self.simage = Instance(image)
        bot.info(self.simage)


def setenv(self, variable, value):
    """set an environment variable for Singularity

    Parameters
    ==========
    variable: the variable to set
    value: the value to set
    """
    os.environ[variable] = value
    os.putenv(variable, value)
    bot.debug("%s set to %s" % (variable, value))


def get_filename(self, image, ext="sif", pwd=True):
    """return an image filename based on the image uri.

    Parameters
    ==========
    ext: the extension to use
    pwd: derive a filename for the pwd
    """
    if pwd:
        image = os.path.basename(image)
    image = re.sub("^.*://", "", image)
    if not image.endswith(ext):
        image = "%s.%s" % (image, ext)
    return image


def get_uri(self):
    """check if the loaded image object (self.simage) has an associated uri
    return if yes, None if not.
    """
    if hasattr(self, "simage"):
        if self.simage is not None:
            if self.simage.image not in ["", None]:
                # Concatenates the <uri>://<image>
                return str(self.simage)
