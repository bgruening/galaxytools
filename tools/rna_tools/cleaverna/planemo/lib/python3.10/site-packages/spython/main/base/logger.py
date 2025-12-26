# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


import os

from spython.logger import decodeUtf8String


def init_level(self, quiet=False):
    """set the logging level based on the environment

    Parameters
    ==========
    quiet: boolean if True, set to quiet. Gets overridden by environment
           setting, and only exists to define default

    """

    if os.environ.get("MESSAGELEVEL") == "QUIET":
        quiet = True

    self.quiet = quiet


def println(self, output, quiet=False):
    """print will print the output, given that quiet is not True. This
    function also serves to convert output in bytes to utf-8

    Parameters
    ==========
    output: the string to print
    quiet: a runtime variable to over-ride the default.

    """
    if not self.quiet and not quiet:
        print(decodeUtf8String(output))
