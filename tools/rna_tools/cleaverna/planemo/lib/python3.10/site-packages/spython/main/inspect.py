# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import json as jsonp

from spython.logger import bot
from spython.utils import check_install, run_command


def inspect(
    self, image=None, json=True, app=None, quiet=True, singularity_options=None
):
    """inspect will show labels, defile, runscript, and tests for an image

    Parameters
    ==========
    image: path of image to inspect
    json: print json instead of raw text (default True)
    quiet: Don't print result to the screen (default True)
    app: if defined, return help in context of an app
    singularity_options: a list of options to provide to the singularity client

    """
    check_install()

    # No image provided, default to use the client's loaded image
    if not image:
        image = self._get_uri()

    # If there still isn't an image, exit on error
    if not image:
        bot.exit("Please provide an image to inspect.")

    cmd = self._init_command("inspect", singularity_options)
    if app:
        cmd = cmd + ["--app", app]

    options = ["e", "d", "l", "r", "H", "t"]
    for x in options:
        cmd.append("-%s" % x)

    if json:
        cmd.append("--json")

    cmd.append(image)

    # Does the user want to see the command printed?
    if not (quiet or self.quiet):
        bot.info(" ".join(cmd))

    result = run_command(cmd, quiet=quiet)

    if result["return_code"] == 0:
        result = jsonp.loads(result["message"][0])

        # Unify output to singularity 3 format
        if "data" in result:
            result = result["data"]

        # Fix up labels
        result = parse_labels(result)

        if not quiet:
            print(jsonp.dumps(result, indent=4))

    return result


def parse_labels(result):
    """fix up the labels, meaning parse to json if needed, and return
    original updated object

    Parameters
    ==========
    result: the json object to parse from inspect
    """

    labels = result["attributes"].get("labels") or {}
    try:
        labels = jsonp.loads(labels)
    except Exception:
        pass

    result["attributes"]["labels"] = labels

    return result
