# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import json

from spython.logger import bot
from spython.utils import stream_command


def run(
    self,
    image=None,
    args=None,
    app=None,
    sudo=False,
    writable=False,
    contain=False,
    bind=None,
    stream=False,
    nv=False,
    options=None,
    singularity_options=None,
    return_result=False,
    quiet=False,
    background=False,
    stream_type="stdout",
):
    """
    run will run the container, with or withour arguments (which
    should be provided in a list)

    Parameters
    ==========
    image: full path to singularity image
    args: args to include with the run
    app: if not None, execute a command in context of an app
    writable: This option makes the file system accessible as read/write
    options: an optional list of options to provide to run.
    singularity_options: a list of options to provide to the singularity client
    contain: This option disables the automatic sharing of writable
             filesystems on your host
    bind: list or single string of bind paths.
          This option allows you to map directories on your host system to
          directories within your container using bind mounts
    stream: if True, return <generator> for the user to run
    nv: if True, load Nvidia Drivers in runtime (default False)
    return_result: if True, return entire json object with return code
         and message result (default is False)
    quiet: print the command to the user
    stream_type: Sets which output stream from the singularity command should be return. Values are 'stdout', 'stderr', 'both'.
    """
    from spython.utils import check_install

    check_install()

    cmd = self._init_command("run", singularity_options)

    # Does the user want to see the command printed?
    quiet = quiet or self.quiet

    # nv option leverages any GPU cards
    if nv:
        cmd += ["--nv"]

    # No image provided, default to use the client's loaded image
    if image is None:
        image = self._get_uri()

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

    # Does the user want writable?
    if writable:
        cmd.append("--writable")

    # Add options
    if options is not None:
        cmd = cmd + options

    cmd = cmd + [image]

    if args is not None:
        if not isinstance(args, list):
            args = args.split(" ")
        cmd = cmd + args

    if not quiet:
        bot.info(" ".join(cmd))

    if background:
        return self._run_command(cmd, sudo=sudo, background=True)

    elif not stream:
        result = self._run_command(cmd, sudo=sudo, return_result=return_result)
    else:
        return stream_command(cmd, sudo=sudo, output_type=stream_type)

    # If the user wants the raw result object
    if return_result:
        return result

    # Otherwise, we parse the result if it was successful
    if result:
        result = result.strip("\n")

        try:
            result = json.loads(result)
        except Exception:
            pass
        return result
