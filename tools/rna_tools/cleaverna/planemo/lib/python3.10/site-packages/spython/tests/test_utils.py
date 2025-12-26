#!/usr/bin/python

# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.

import os

import pytest

from spython.utils import ScopedEnvVar


def test_write_read_files(tmp_path):
    """
    test_write_read_files will test the functions write_file and read_file
    """
    print("Testing utils.write_file...")
    from spython.utils import write_file

    tmpfile = str(tmp_path / "written_file.txt")
    assert not os.path.exists(tmpfile)
    write_file(tmpfile, "hello!")
    assert os.path.exists(tmpfile)

    print("Testing utils.read_file...")
    from spython.utils import read_file

    content = read_file(tmpfile)[0]
    assert content == "hello!"


def test_write_bad_json(tmp_path):
    from spython.utils import write_json

    bad_json = {"Wakkawakkawakka'}": [{True}, "2", 3]}
    tmpfile = str(tmp_path / "json_file.txt")
    assert not os.path.exists(tmpfile)
    with pytest.raises(TypeError):
        write_json(bad_json, tmpfile)


def test_write_json(tmp_path):
    import json

    from spython.utils import write_json

    good_json = {"Wakkawakkawakka": [True, "2", 3]}
    tmpfile = str(tmp_path / "good_json_file.txt")
    assert not os.path.exists(tmpfile)
    write_json(good_json, tmpfile)
    with open(tmpfile, "r") as f:
        content = json.loads(f.read())
    assert isinstance(content, dict)
    assert "Wakkawakkawakka" in content


def test_check_install():
    """check install is used to check if a particular software is installed.
    If no command is provided, singularity is assumed to be the test case"""
    print("Testing utils.check_install")
    from spython.utils import check_install

    is_installed = check_install()
    assert is_installed
    is_not_installed = check_install("fakesoftwarename")
    assert not is_not_installed


def test_check_get_singularity_version():
    """check that the singularity version is found to be that installed"""
    from spython.utils import get_singularity_version

    version = get_singularity_version()
    assert version != ""
    with ScopedEnvVar("SPYTHON_SINGULARITY_VERSION", "3.0"):
        version = get_singularity_version()
    assert version == "3.0"


def test_get_installdir():
    """get install directory should return the base of where singularity
    is installed
    """
    print("Testing utils.get_installdir")
    from spython.utils import get_installdir

    whereami = get_installdir()
    print(whereami)
    assert whereami.endswith("spython")


def test_split_uri():
    from spython.utils import split_uri

    protocol, image = split_uri("docker://ubuntu")
    assert protocol == "docker"
    assert image == "ubuntu"

    protocol, image = split_uri("http://image/path/with/slash/")
    assert protocol == "http"
    assert image == "image/path/with/slash"

    protocol, image = split_uri("no/proto/")
    assert protocol == ""
    assert image == "no/proto"


def test_remove_uri():
    print("Testing utils.remove_uri")
    from spython.utils import remove_uri

    assert remove_uri("docker://ubuntu") == "ubuntu"
    assert (
        remove_uri("shub://vanessa/singularity-images") == "vanessa/singularity-images"
    )
    assert remove_uri("library://library/default/alpine") == "library/default/alpine"
    assert remove_uri("vanessa/singularity-images") == "vanessa/singularity-images"


def test_decode():
    from spython.logger import decodeUtf8String

    out = decodeUtf8String(str("Hello"))
    assert isinstance(out, str)
    assert out == "Hello"
    out = decodeUtf8String(bytes(b"Hello"))
    assert isinstance(out, str)
    assert out == "Hello"


def test_ScopedEnvVar():
    assert "FOO" not in os.environ
    with ScopedEnvVar("FOO", "bar") as e:
        assert e.name == "FOO"
        assert e.value == "bar"
        assert os.environ["FOO"] == "bar"
        with ScopedEnvVar("FOO", "baz"):
            assert os.environ["FOO"] == "baz"
        assert os.environ["FOO"] == "bar"
        # None removes it
        with ScopedEnvVar("FOO", None):
            assert "FOO" not in os.environ
        # But empty string is allowed
        with ScopedEnvVar("FOO", ""):
            assert os.environ["FOO"] == ""
        assert os.environ["FOO"] == "bar"
    assert "FOO" not in os.environ
    # Unset a non-existing variable
    with ScopedEnvVar("FOO", None):
        assert "FOO" not in os.environ
    assert "FOO" not in os.environ
