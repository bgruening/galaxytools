# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import shutil

from spython.logger import bot
from spython.utils import stream_command


def execute(
    self,
    image=None,
    command=None,
    app=None,
    writable=False,
    contain=False,
    bind=None,
    stream=False,
    nv=False,
    return_result=False,
    options=None,
    singularity_options=None,
    sudo=False,
    sudo_options=None,
    quiet=True,
    environ=None,
    stream_type="stdout",
):
    """execute: send a command to a container

    Parameters
    ==========

    image: full path to singularity image
    command: command to send to container
    app: if not None, execute a command in context of an app
    writable: This option makes the file system accessible as read/write
    contain: This option disables the automatic sharing of writable
             filesystems on your host
    options: an optional list of options to provide to execute.
    singularity_options: a list of options to provide to the singularity client
    bind: list or single string of bind paths.
         This option allows you to map directories on your host system to
         directories within your container using bind mounts
    nv: if True, load Nvidia Drivers in runtime (default False)
    return_result: if True, return entire json object with return code
                   and message result not (default)
    quiet: Do not print verbose output.
    environ: extra environment to add.
    stream_type: Sets which output stream from the singularity command should be return. Values are 'stdout', 'stderr', 'both'.
    """
    from spython.utils import check_install

    check_install()

    cmd = self._init_command("exec", singularity_options)

    # nv option leverages any GPU cards
    if nv:
        cmd += ["--nv"]

    # If the image is given as a list, it's probably the command
    if isinstance(image, list):
        command = image
        image = None

    if command is not None:
        # No image provided, default to use the client's loaded image
        if image is None:
            image = self._get_uri()
            self.quiet = True

        # If an instance is provided, grab it's name
        if isinstance(image, self.instance):
            image = image.get_uri()

        # If image is still None, not defined by user or previously with client
        if image is None:
            bot.exit("Please load or provide an image.")

        # Does the user want to use bind paths option?
        if bind is not None:
            cmd += self._generate_bind_list(bind)

        # Does the user want to run an app?
        if app is not None:
            cmd = cmd + ["--app", app]

        if writable:
            cmd.append("--writable")

        # Add additional options
        if options is not None:
            cmd = cmd + options

        if not isinstance(command, list):
            command = command.split(" ")

        cmd = cmd + [image] + command

        # Does the user want to see the command printed?
        if not (quiet or self.quiet):
            bot.info(" ".join(cmd))

        if not stream:
            return self._run_command(
                cmd,
                sudo=sudo,
                sudo_options=sudo_options,
                return_result=return_result,
                quiet=quiet,
                environ=environ,
            )
        return stream_command(
            cmd, sudo=sudo, sudo_options=sudo_options, output_type=stream_type
        )

    bot.exit("Please include a command (list) to execute.")


def shell(
    self,
    image,
    app=None,
    writable=False,
    contain=False,
    bind=None,
    nv=False,
    options=None,
    singularity_options=None,
    sudo=False,
    quiet=True,
):
    """shell into a container. A user is advised to use singularity to do
    this directly, however this function is useful for supporting tools.

    Parameters
    ==========

    image: full path to singularity image
    app: if not None, execute a shell in context of an app
    writable: This option makes the file system accessible as read/write
    contain: This option disables the automatic sharing of writable
             filesystems on your host
    options: an optional list of options to provide to shell.
    singularity_options: a list of options to provide to the singularity client
    bind: list or single string of bind paths.
         This option allows you to map directories on your host system to
         directories within your container using bind mounts
    nv: if True, load Nvidia Drivers in runtime (default False)
    """
    from spython.utils import check_install

    check_install()

    cmd = self._init_command("shell", singularity_options)

    # nv option leverages any GPU cards
    if nv:
        cmd += ["--nv"]

    # Does the user want to use bind paths option?
    if bind is not None:
        cmd += self._generate_bind_list(bind)

    # Does the user want to run an app?
    if app is not None:
        cmd = cmd + ["--app", app]

    # Add additional options
    if options is not None:
        cmd = cmd + options

    if writable:
        cmd.append("--writable")

    # Finally, add the image or uri
    cmd.append(image)
    singularity = shutil.which("singularity")
    if not singularity:
        raise ValueError("Cannot find singularity executable.")

    # Does the user want to see the command printed?
    if not (quiet or self.quiet):
        bot.info(" ".join(cmd))

    if writable or sudo:
        os.execvp("sudo", ["sudo"] + cmd)

    else:
        os.execvp(singularity, cmd)
