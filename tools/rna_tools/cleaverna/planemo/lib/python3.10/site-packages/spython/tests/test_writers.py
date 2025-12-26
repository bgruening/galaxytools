#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os

from spython.main.parse.writers import DockerWriter, SingularityWriter


def test_writers():
    from spython.main.parse.writers import get_writer

    writer = get_writer("docker")
    assert writer == DockerWriter

    writer = get_writer("Dockerfile")
    assert writer == DockerWriter

    writer = get_writer("Singularity")
    assert writer == SingularityWriter


def test_docker_writer(test_data):
    from spython.main.parse.parsers import DockerParser

    dockerfile = os.path.join(test_data["root"], "Dockerfile")
    parser = DockerParser(dockerfile)
    writer = DockerWriter(parser.recipe)

    assert str(writer) == "[spython-writer][docker]"
    print(writer.convert())


def test_singularity_writer(test_data):
    from spython.main.parse.parsers import SingularityParser

    recipe = os.path.join(test_data["root"], "Singularity")
    parser = SingularityParser(recipe)
    writer = SingularityWriter(parser.recipe)

    assert str(writer) == "[spython-writer][singularity]"
    print(writer.convert())
