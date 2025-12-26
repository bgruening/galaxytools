# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import re

from spython.logger import bot
from spython.utils import ScopedEnvVar, stream_command


def pull(
    self,
    image=None,
    name=None,
    pull_folder="",
    ext="sif",
    force=False,
    capture=False,
    stream=False,
    quiet=False,
    singularity_options=None,
):
    """pull will pull a singularity hub or Docker image

    Parameters
    ==========
    image: the complete image uri. If not provided, the client loaded is used
    singularity_options: a list of options to provide to the singularity client
    pull_folder: if not defined, pulls to $PWD (''). If defined, pulls to
                 user specified location instead.

    Docker and Singularity Hub Naming
    ---------------------------------
    name: a custom name to use, to override default
    ext: if no name specified, the default extension to use.

    """
    from spython.utils import check_install

    check_install()

    cmd = self._init_command("pull", singularity_options)

    # Quiet is honored if set by the client, or user
    quiet = quiet or self.quiet

    # No image provided, default to use the client's loaded image
    if image is None:
        image = self._get_uri()

    # If it's still None, no go!
    if image is None:
        bot.exit("You must provide an image uri, or use client.load() first.")

    # Singularity Only supports shub, docker and library pull
    if not re.search("^(shub|docker|library|https|oras)://", image):
        bot.exit("pull only valid for docker, oras, https, shub and library.")

    # If we still don't have a custom name, base off of image uri.
    if name is None:
        name = self._get_filename(image, ext)

    if pull_folder:
        final_image = os.path.join(pull_folder, os.path.basename(name))

        # Regression Singularity 3.* onward, PULLFOLDER not honored
        # https://github.com/sylabs/singularity/issues/2788
        name = final_image
        pull_folder = None  # Don't use pull_folder
    else:
        final_image = name

    cmd = cmd + ["--name", name]

    if force:
        cmd = cmd + ["--force"]

    cmd.append(image)

    if not quiet:
        bot.info(" ".join(cmd))

    with ScopedEnvVar("SINGULARITY_PULLFOLDER", pull_folder):
        # Option 1: Streaming we just run to show user
        if not stream:
            self._run_command(cmd, capture=capture, quiet=quiet)

        # Option 3: A custom name we can predict (not commit/hash) and can also show
        else:
            # As of Singularity 3.x (at least 3.8) output goes to stderr
            return final_image, stream_command(cmd, sudo=False, output_type="stderr")

    if os.path.exists(final_image) and not quiet:
        bot.info(final_image)
    return final_image
