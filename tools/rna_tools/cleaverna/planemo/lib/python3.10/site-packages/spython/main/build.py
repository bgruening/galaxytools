# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import re

from spython.logger import bot
from spython.utils import stream_command


def build(
    self,
    recipe=None,
    image=None,
    isolated=False,
    sandbox=False,
    writable=False,
    build_folder=None,
    robot_name=False,
    ext="sif",
    sudo=True,
    stream=False,
    force=False,
    options=None,
    quiet=False,
    return_result=False,
    sudo_options=None,
    singularity_options=None,
):
    """build a singularity image, optionally for an isolated build
    (requires sudo). If you specify to stream, expect the image name
    and an iterator to be returned.

    image, builder = Client.build(...)

    Parameters
    ==========

    recipe: the path to the recipe file (or source to build from). If not
               defined, we look for "Singularity" file in $PWD
    image: the image to build (if None, will use arbitrary name
    isolated: if True, run build with --isolated flag
    sandbox: if True, create a writable sandbox
    writable: if True, use writable ext3 (sandbox takes preference)
    build_folder: where the container should be built.
    ext: the image extension to use.
    robot_name: boolean, default False. if you don't give your image a
                name (with "image") then a fun robot name will be generated
                instead. Highly recommended :)
    sudo: give sudo to the command (or not) default is True for build
    sudo_options: options to pass to sudo (e.g. --preserve-env=SINGULARITY_CACHEDIR,SINGULARITY_TMPDIR)
    options: for all other options, specify them in this list.
    singularity_options: a list of options to provide to the singularity client
    quiet: quiet verbose printing from the client.
    return_result: if True, return complete error code / message dictionary
    """
    from spython.utils import check_install

    check_install()

    cmd = self._init_command("build", singularity_options)

    # If no extra options
    options = options or []

    # Force the build if the image / sandbox exists
    if force:
        cmd.append("--force")

    # No image provided, default to use the client's loaded image
    if recipe is None:
        recipe = self._get_uri()

    # If it's still None, try default build recipe
    if recipe is None:
        recipe = "Singularity"

        if not os.path.exists(recipe):
            bot.exit("Cannot find %s, exiting." % image)

    if image is None:
        if re.search("(docker|shub|library)://", recipe) and not robot_name:
            image = self._get_filename(recipe, ext)
        else:
            image = "%s.%s" % (self.RobotNamer.generate(), ext)

    # Does the user want a custom build folder?
    if build_folder is not None:
        if not os.path.exists(build_folder):
            bot.exit("%s does not exist!" % build_folder)
        image = os.path.join(build_folder, image)

    # The user wants to run an isolated build
    if isolated:
        cmd.append("--isolated")

    if sandbox:
        cmd.append("--sandbox")

    elif writable:
        cmd.append("--writable")

    cmd = cmd + options + [image, recipe]

    # Does the user want to see the command printed?
    if not (quiet or self.quiet):
        bot.info(" ".join(cmd))

    if not stream:
        self._run_command(
            cmd,
            sudo=sudo,
            sudo_options=sudo_options,
            quiet=quiet,
            return_result=return_result,
            capture=False,
        )

    else:
        # Here we return the expected image, and an iterator!
        # The caller must iterate over
        return image, stream_command(cmd, sudo=sudo, sudo_options=sudo_options)

    if os.path.exists(image):
        return image
