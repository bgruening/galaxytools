import os
from glob import glob

import pytest

from spython.main import Client
from spython.utils import get_installdir


@pytest.fixture
def installdir():
    return get_installdir()


@pytest.fixture
def test_data(installdir):  # pylint: disable=redefined-outer-name
    root = os.path.join(installdir, "tests", "testdata")
    dockerFiles = glob(os.path.join(root, "docker2singularity", "*.docker"))
    singularityFiles = glob(os.path.join(root, "singularity2docker", "*.def"))
    return {
        "root": root,
        "d2s": [(file, os.path.splitext(file)[0] + ".def") for file in dockerFiles],
        "s2d": [
            (file, os.path.splitext(file)[0] + ".docker") for file in singularityFiles
        ],
    }


@pytest.fixture(scope="session")
def oras_container(tmp_path_factory):
    folder = tmp_path_factory.mktemp("oras-img")
    return folder, Client.pull(
        "oras://ghcr.io/singularityhub/github-ci:latest", pull_folder=str(folder)
    )


@pytest.fixture(scope="session")
def docker_container(tmp_path_factory):
    folder = tmp_path_factory.mktemp("docker-img")
    return folder, Client.pull("docker://busybox:1.30.1", pull_folder=str(folder))
