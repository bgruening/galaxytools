# Copyright (C) 2019-2022 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def mount(self, image, sudo=None):
    """create an OCI bundle from SIF image

    Parameters
    ==========
    image: the container (sif) to mount
    """
    return self._state_command(image, command="mount", sudo=sudo)


def umount(self, image, sudo=None):
    """delete an OCI bundle created from SIF image

    Parameters
    ==========
    image: the container (sif) to mount
    """
    return self._state_command(image, command="umount", sudo=sudo)
