#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def test_recipe_base():
    from spython.main.parse.recipe import Recipe

    recipe = Recipe()
    assert str(recipe) == "[spython-recipe]"

    attributes = [
        "cmd",
        "comments",
        "entrypoint",
        "environ",
        "files",
        "install",
        "labels",
        "ports",
        "test",
        "volumes",
        "workdir",
    ]

    for att in attributes:
        assert hasattr(recipe, att)

    print("Checking that empty recipe returns empty")
    result = recipe.json()
    assert not result

    print("Checking that non-empty recipe returns values")
    recipe.cmd = ["echo", "hello"]
    recipe.entrypoint = "/bin/bash"
    recipe.comments = ["This recipe is great", "Yes it is!"]
    recipe.environ = ["PANCAKES=WITHSYRUP"]
    recipe.files = [["one", "two"]]
    recipe.test = ["true"]
    recipe.install = ["apt-get update"]
    recipe.labels = ["Maintainer vanessasaur"]
    recipe.ports = ["3031"]
    recipe.volumes = ["/data"]
    recipe.workdir = "/code"

    result = recipe.json()
    for att in attributes:
        assert att in result
