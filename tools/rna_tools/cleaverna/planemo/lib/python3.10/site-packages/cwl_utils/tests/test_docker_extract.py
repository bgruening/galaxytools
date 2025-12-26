# SPDX-License-Identifier: Apache-2.0
"""Tests for cwl-docker-extract."""
from pathlib import Path

import pytest

from cwl_utils.docker_extract import arg_parser, run

from .util import get_data, needs_docker, needs_podman, needs_singularity


@pytest.mark.parametrize(
    ("target", "engine"),
    [
        pytest.param("testdata/md5sum.cwl", "docker", marks=needs_docker),
        pytest.param("testdata/md5sum_v11.cwl", "docker", marks=needs_docker),
        pytest.param("testdata/md5sum.cwl", "podman", marks=needs_podman),
        pytest.param("testdata/md5sum_v11.cwl", "podman", marks=needs_podman),
        pytest.param("testdata/md5sum.cwl", "singularity", marks=needs_singularity),
        pytest.param("testdata/md5sum_v11.cwl", "singularity", marks=needs_singularity),
    ],
)
def test_container_extraction(target: str, engine: str, tmp_path: Path) -> None:
    """Test container extraction tool."""
    args = ["--dir", str(tmp_path), get_data(target), "--container-engine", engine]
    if engine == "singularity":
        args.append("--singularity")
    reqs = run(arg_parser().parse_args(args))
    assert len(reqs) == 1
    assert len(list(tmp_path.iterdir())) == 1


@pytest.mark.parametrize(
    ("engine"),
    [
        pytest.param("docker", marks=needs_docker),
        pytest.param("podman", marks=needs_podman),
        pytest.param("singularity", marks=needs_singularity),
    ],
)
def test_container_extraction_force(engine: str, tmp_path: Path) -> None:
    """Test force pull container extraction."""
    args = [
        "--dir",
        str(tmp_path),
        get_data("testdata/md5sum.cwl"),
        "--container-engine",
        engine,
    ]
    if engine == "singularity":
        args.append("--singularity")
    reqs = run(arg_parser().parse_args(args))
    assert len(reqs) == 1
    assert len(list(tmp_path.iterdir())) == 1
    args = [
        "--dir",
        str(tmp_path),
        get_data("testdata/md5sum.cwl"),
        "--container-engine",
        engine,
        "--force-download",
    ]
    if engine == "singularity":
        args.append("--singularity")
    reqs = run(arg_parser().parse_args(args))
    assert len(reqs) == 1
    assert len(list(tmp_path.iterdir())) == 1


@pytest.mark.parametrize(
    ("engine"),
    [
        pytest.param("docker", marks=needs_docker),
        pytest.param("podman", marks=needs_podman),
        pytest.param("singularity", marks=needs_singularity),
    ],
)
def test_container_extraction_no_dockerPull(
    engine: str, tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    """Test container extraction tool when dockerPull is missing."""
    args = [
        "--dir",
        str(tmp_path),
        get_data("testdata/debian_image_id.cwl"),
        "--container-engine",
        engine,
    ]
    if engine == "singularity":
        args.append("--singularity")
    reqs = run(arg_parser().parse_args(args))
    assert len(reqs) == 1
    assert len(list(tmp_path.iterdir())) == 0
    captured = capsys.readouterr()
    assert (
        captured.err
        == """Unable to save image from due to lack of 'dockerPull':
class: DockerRequirement
dockerImageId: 'debian:stable-slim.img'
"""
    )


@pytest.mark.parametrize(
    ("engine"),
    [
        pytest.param("docker", marks=needs_docker),
        pytest.param("podman", marks=needs_podman),
        pytest.param("singularity", marks=needs_singularity),
    ],
)
def test_container_extraction_embedded_step(engine: str, tmp_path: Path) -> None:
    """Test container extraction tool."""
    args = [
        "--dir",
        str(tmp_path),
        get_data("testdata/workflows/count-lines16-wf.cwl"),
        "--container-engine",
        engine,
    ]
    if engine == "singularity":
        args.append("--singularity")
    reqs = run(arg_parser().parse_args(args))
    assert len(reqs) == 1
    assert len(list(tmp_path.iterdir())) == 1
