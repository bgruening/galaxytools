# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import hashlib
import os

from spython.logger import bot
from spython.utils import split_uri


class ImageBase:
    def __str__(self):
        protocol = getattr(self, "protocol", None)
        if protocol:
            return "%s://%s" % (protocol, self.image)
        return self.image

    def __repr__(self):
        return self.__str__()

    def parse_image_name(self, image):
        """
        simply split the uri from the image. Singularity handles
        parsing of registry, namespace, image.

        Parameters
        ==========
        image: the complete image uri to load (e.g., docker://ubuntu)

        """
        self._image = image
        self.protocol, self.image = split_uri(image)


class Image(ImageBase):
    def __init__(self, image=None):
        """An image here is an image file or a record.
        The user can choose to load the image when starting the client, or
        update the main client with an image. The image object is kept
        with the main client to make running additional commands easier.

        Parameters
        ==========
        image: the image uri to parse (required)

        """
        super(Image, self).__init__()
        self.parse_image_name(image)

    def get_hash(self, image=None):
        """return an md5 hash of the file based on a criteria level. This
        is intended to give the file a reasonable version. This only is
        useful for actual image files.

        Parameters
        ==========
        image: the image path to get hash for (first priority). Second
               priority is image path saved with image object, if exists.

        """
        hasher = hashlib.md5()
        image = image or self.image

        if os.path.exists(image):
            with open(image, "rb") as f:
                for chunk in iter(lambda: f.read(4096), b""):
                    hasher.update(chunk)
                return hasher.hexdigest()

        bot.warning("%s does not exist." % image)
