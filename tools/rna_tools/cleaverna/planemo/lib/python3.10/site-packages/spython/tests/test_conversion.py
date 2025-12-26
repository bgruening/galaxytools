#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os
from glob import glob


def read_file(file):
    with open(file) as fd:
        content = fd.read().strip("\n")
    return content


def test_other_recipe_exists(test_data):
    # Have any example
    assert test_data["d2s"]
    assert test_data["s2d"]

    for _, outFile in test_data["d2s"] + test_data["s2d"]:
        assert os.path.exists(outFile), outFile + " is missing"

    dockerfiles = glob(
        os.path.join(os.path.dirname(test_data["s2d"][0][0]), "*.docker")
    )
    singularityfiles = glob(
        os.path.join(os.path.dirname(test_data["d2s"][0][0]), "*.def")
    )
    for file in dockerfiles:
        assert file in [out for _, out in test_data["s2d"]]
    for file in singularityfiles:
        assert file in [out for _, out in test_data["d2s"]]


def test_docker2singularity(test_data, tmp_path):
    from spython.main.parse.parsers import DockerParser
    from spython.main.parse.writers import SingularityWriter

    for dockerfile, recipe in test_data["d2s"]:
        parser = DockerParser(dockerfile)
        writer = SingularityWriter(parser.recipe)
        result = read_file(recipe).strip()
        assert writer.convert().replace("\n", "") == result.replace("\n", "")


def test_singularity2docker(test_data, tmp_path):
    print("Testing spython conversion from singularity2docker")
    from spython.main.parse.parsers import SingularityParser
    from spython.main.parse.writers import DockerWriter

    for recipe, dockerfile in test_data["s2d"]:
        parser = SingularityParser(recipe)
        writer = DockerWriter(parser.recipe)
        assert writer.convert() == read_file(dockerfile)
