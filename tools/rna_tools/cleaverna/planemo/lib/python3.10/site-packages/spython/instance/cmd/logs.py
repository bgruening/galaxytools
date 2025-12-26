# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import os
import platform

from spython.logger import bot
from spython.utils import get_userhome, get_username


def error_logs(self, print_logs=False):
    """For Singularity 3.5 and later, we are able to programatically
    derive the name of the log. In this case, return the content
    to the user. See
    https://github.com/sylabs/singularity/issues/1115#issuecomment-560457918
    for when this was added.

    Parameters
    ==========
    print_logs: boolean to indicate to print to the screen along with
                return (defaults to False to just return log string)
    """
    return self._logs(print_logs, "err")


def output_logs(self, print_logs=False):
    """Get output logs for the user, if they exist.

    Parameters
    ==========
    print_logs: boolean to indicate to print to the screen along with
                return (defaults to False to just return log string)
    """
    return self._logs(print_logs, "out")


def _logs(self, print_logs=False, ext="out"):
    """A shared function to print log files. The only differing element is
    the extension (err or out)
    """
    from spython.utils import check_install

    check_install()

    # Formulate the path of the logs
    hostname = platform.node()
    logpath = os.path.join(
        get_userhome(),
        ".singularity",
        "instances",
        "logs",
        hostname,
        get_username(),
        "%s.%s" % (self.name, ext),
    )

    if os.path.exists(logpath):
        with open(logpath, "r") as filey:
            logs = filey.read()
        if print_logs is True:
            print(logs)
    else:
        bot.warning("No log files have been produced.")
    return logs
