#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

from spython.image import Image


def test_image():
    image = Image("docker://ubuntu")
    assert str(image) == "docker://ubuntu"
    assert image.protocol == "docker"
    assert image.image == "ubuntu"
