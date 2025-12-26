# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os


def setEnvVar(name, value):
    """Set or unset an environment variable

    name -- Name of the variable to set
    value -- Value to use or None to clear
    """
    if value is None:
        if name in os.environ:
            del os.environ[name]
    else:
        os.environ[name] = value


class ScopedEnvVar:
    """Temporarily change an environment variable

    Usage:
        with ScopedEnvVar("FOO", "bar"):
            print(os.environ["FOO"]) # "bar"
        print(os.environ["FOO"]) # <oldvalue>
    """

    def __init__(self, name, value):
        """Create the scoped environment variable object

        name -- Name of the variable to set
        value -- Value to use or None to clear
        """
        self.name = name
        self.value = value
        self.oldValue = None

    def __enter__(self):
        self.oldValue = os.environ.get(self.name)
        setEnvVar(self.name, self.value)
        return self

    def __exit__(self, ex_type, ex_value, traceback):
        setEnvVar(self.name, self.oldValue)
