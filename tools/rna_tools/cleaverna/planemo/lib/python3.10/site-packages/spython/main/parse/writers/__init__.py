# Copyright (C) 2017-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


from .docker import DockerWriter
from .singularity import SingularityWriter


def get_writer(name):
    """get_writer is a simple helper function to return a writer based on it's
    name, if it exists. If there is no writer defined, we return None.

    Parameters
    ==========
    name: the name of the writer to return.
    """
    name = name.lower()
    writers = {
        "docker": DockerWriter,
        "singularity": SingularityWriter,
        "dockerfile": DockerWriter,
    }

    if name in writers:
        return writers[name]
