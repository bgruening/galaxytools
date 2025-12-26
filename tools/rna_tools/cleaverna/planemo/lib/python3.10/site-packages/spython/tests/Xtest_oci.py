#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
import shutil

import pytest

from spython.main import Client
from spython.main.base.generate import RobotNamer
from spython.utils import get_installdir


@pytest.fixture
def sandbox(tmp_path):
    image = Client.build(
        "docker://busybox:1.30.1",
        image=str(tmp_path / "sandbox"),
        sandbox=True,
        sudo=False,
    )

    assert os.path.exists(image)

    config = os.path.join(get_installdir(), "oci", "config.json")
    shutil.copyfile(config, os.path.join(image, "config.json"))
    return image


def test_oci_image():
    image = Client.oci.OciImage("oci://imagename")
    assert image.get_uri() == "[singularity-python-oci:oci://imagename]"


def test_oci(sandbox):  # pylint: disable=redefined-outer-name
    image = sandbox
    container_id = RobotNamer().generate()

    # A non existing process should not have a state
    print("...Case 1. Check status of non-existing bundle.")
    state = Client.oci.state("mycontainer")
    assert state is None

    # This will use sudo
    print("...Case 2: Create OCI image from bundle")
    result = Client.oci.create(bundle=image, container_id=container_id)

    print(result)
    assert result["status"] == "created"

    print("...Case 3. Execute command to non running bundle.")
    result = Client.oci.execute(
        container_id=container_id, sudo=True, command=["ls", "/"]
    )

    print(result)
    assert "bin" in result

    print("...Case 4. Start container return value 0.")
    state = Client.oci.start(container_id, sudo=True)
    assert state == 0

    print("...Case 5. Execute command to running bundle.")
    result = Client.oci.execute(
        container_id=container_id, sudo=True, command=["ls", "/"]
    )

    print(result)
    assert "bin" in result

    print("...Case 6. Check status of existing bundle.")
    state = Client.oci.state(container_id, sudo=True)
    assert state["status"] == "running"

    print("...Case 7. Pause running container return value 0.")
    state = Client.oci.pause(container_id, sudo=True)
    assert state == 0

    # State was still reported as running
    print("...check status of paused bundle.")
    state = Client.oci.state(container_id, sudo=True)
    assert state["status"] == "paused"

    print("...Case 8. Resume paused container return value 0.")
    state = Client.oci.resume(container_id, sudo=True)
    assert state == 0

    print("...check status of resumed bundle.")
    state = Client.oci.state(container_id, sudo=True)
    assert state["status"] == "running"

    print("...Case 9. Kill container.")
    state = Client.oci.kill(container_id, sudo=True)
    assert state == 0

    # Clean up the image (should still use sudo)
    # Bug in singularity that kill doesn't kill completely - this returns
    # 255. When testsupdated to 3.1.* add signal=K to run
    result = Client.oci.delete(container_id, sudo=True)
    assert result in [0, 255]
