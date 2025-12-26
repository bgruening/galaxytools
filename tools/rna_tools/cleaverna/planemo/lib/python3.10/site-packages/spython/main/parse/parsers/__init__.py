# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

from .docker import DockerParser
from .singularity import SingularityParser


def get_parser(name):
    """get_parser is a simple helper function to return a parser based on it's
    name, if it exists. If there is no writer defined, we return None.

    Parameters
    ==========
    name: the name of the parser to return.
    """
    name = name.lower()
    parsers = {
        "docker": DockerParser,
        "singularity": SingularityParser,
        "dockerfile": DockerParser,
    }

    if name in parsers:
        return parsers[name]
