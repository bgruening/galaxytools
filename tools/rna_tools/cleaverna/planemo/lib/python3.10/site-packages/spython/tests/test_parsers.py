#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os

from spython.main.parse.parsers import DockerParser, SingularityParser


def test_get_parser():
    from spython.main.parse.parsers import get_parser

    parser = get_parser("docker")
    assert parser == DockerParser

    parser = get_parser("Dockerfile")
    assert parser == DockerParser

    parser = get_parser("Singularity")
    assert parser == SingularityParser


def test_docker_parser(test_data):
    dockerfile = os.path.join(test_data["root"], "Dockerfile")
    parser = DockerParser(dockerfile)

    assert str(parser) == "[spython-parser][docker]"
    assert "spython-base" in parser.recipe
    recipe = parser.recipe["spython-base"]

    # Test all fields from recipe
    assert recipe.fromHeader == "python:3.5.1"
    assert recipe.cmd == "/code/run_uwsgi.sh"
    assert recipe.entrypoint is None
    assert recipe.workdir == "/code"
    assert recipe.volumes == []
    assert recipe.ports == ["3031"]
    assert recipe.files[0] == ["requirements.txt", "/tmp/requirements.txt"]
    assert recipe.environ == ["PYTHONUNBUFFERED=1"]
    assert recipe.source == dockerfile


def test_singularity_parser(test_data):
    recipefile = os.path.join(test_data["root"], "Singularity")
    parser = SingularityParser(recipefile)

    assert str(parser) == "[spython-parser][singularity]"
    assert "spython-base" in parser.recipe
    recipe = parser.recipe["spython-base"]

    # Test all fields from recipe
    assert recipe.fromHeader == "continuumio/miniconda3"
    assert recipe.cmd == 'exec /opt/conda/bin/spython "$@"'
    assert recipe.entrypoint is None
    assert recipe.workdir is None
    assert recipe.volumes == []
    assert recipe.files == []
    assert recipe.environ == []
    assert recipe.source == recipefile
