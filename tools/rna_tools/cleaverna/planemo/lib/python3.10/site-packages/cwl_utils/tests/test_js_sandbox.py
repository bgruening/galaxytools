"""Test sandboxjs.py and related code."""

import os
import shutil
import threading
from pathlib import Path
from typing import Any

import pytest

from cwl_utils import expression, sandboxjs

from .util import needs_podman, needs_singularity, needs_udocker

node_versions = [
    ("v0.8.26\n", False),
    ("v0.10.25\n", False),
    ("v0.10.26\n", True),
    ("v4.4.2\n", True),
    ("v7.7.3\n", True),
]


@pytest.mark.parametrize("version,supported", node_versions)
def test_node_version(version: str, supported: bool, mocker: Any) -> None:
    mocked_subprocess = mocker.patch("cwl_utils.sandboxjs.subprocess")
    mocked_subprocess.check_output = mocker.Mock(return_value=version)

    assert sandboxjs.check_js_threshold_version("node") == supported


def test_value_from_two_concatenated_expressions() -> None:
    js_engine = sandboxjs.get_js_engine()
    js_engine.have_node_slim = False  # type: ignore[attr-defined]
    js_engine.localdata = threading.local()  # type: ignore[attr-defined]
    assert (
        expression.do_eval(
            '$("a ")$("string")',
            {},
            [{"class": "InlineJavascriptRequirement"}],
            None,
            None,
            {},
            cwlVersion="v1.0",
        )
        == "a string"
    )


def hide_nodejs(temp_dir: Path) -> str:
    """Generate a new PATH that hides node{js,}."""
    paths: list[str] = os.environ.get("PATH", "").split(":")
    names: list[str] = []
    if "/bin" in paths:
        bin_path = Path("/bin")
        if (
            bin_path.is_symlink()
            and os.readlink(bin_path) == "usr/bin"
            and "/usr/bin" in paths
        ):
            paths.remove("/bin")
    for name in ("nodejs", "node"):
        path = shutil.which(name, path=":".join(paths))
        if path:
            names.append(path)
    for name in names:
        dirname = os.path.dirname(name)
        if dirname in paths:
            paths.remove(dirname)
            new_dir = temp_dir / os.path.basename(dirname)
            new_dir.mkdir()
            for entry in os.listdir(dirname):
                if entry not in ("nodejs", "node"):
                    os.symlink(os.path.join(dirname, entry), new_dir / entry)
            paths.append(str(new_dir))
    return ":".join(paths)


@needs_podman
def test_value_from_two_concatenated_expressions_podman(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Javascript test using podman."""
    new_paths = hide_nodejs(tmp_path)
    with monkeypatch.context() as m:
        m.setenv("PATH", new_paths)
        js_engine = sandboxjs.get_js_engine()
        js_engine.have_node_slim = False  # type: ignore[attr-defined]
        js_engine.localdata = threading.local()  # type: ignore[attr-defined]
        assert (
            expression.do_eval(
                '$("a ")$("string")',
                {},
                [{"class": "InlineJavascriptRequirement"}],
                None,
                None,
                {},
                cwlVersion="v1.0",
                container_engine="podman",
            )
            == "a string"
        )


@needs_udocker
def test_value_from_two_concatenated_expressions_udocker(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Javascript test using udocker."""
    new_paths = hide_nodejs(tmp_path)
    with monkeypatch.context() as m:
        m.setenv("PATH", new_paths)
        js_engine = sandboxjs.get_js_engine()
        js_engine.have_node_slim = False  # type: ignore[attr-defined]
        js_engine.localdata = threading.local()  # type: ignore[attr-defined]
        assert (
            expression.do_eval(
                '$("a ")$("string")',
                {},
                [{"class": "InlineJavascriptRequirement"}],
                None,
                None,
                {},
                cwlVersion="v1.0",
                container_engine="udocker",
            )
            == "a string"
        )


@needs_singularity
def test_value_from_two_concatenated_expressions_singularity(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    """Javascript test using Singularity."""
    new_paths = hide_nodejs(tmp_path)
    with monkeypatch.context() as m:
        m.setenv("PATH", new_paths)
        js_engine = sandboxjs.get_js_engine()
        js_engine.have_node_slim = False  # type: ignore[attr-defined]
        js_engine.localdata = threading.local()  # type: ignore[attr-defined]
        assert (
            expression.do_eval(
                '$("a ")$("string")',
                {},
                [{"class": "InlineJavascriptRequirement"}],
                None,
                None,
                {},
                cwlVersion="v1.0",
                container_engine="singularity",
            )
            == "a string"
        )


@needs_singularity
def test_singularity_cache(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    """Affirm that CWL_SINGULARIT_CACHE is respected."""
    bin_path = tmp_path / "no_nodejs"
    bin_path.mkdir()
    new_paths = hide_nodejs(bin_path)
    cache_path = tmp_path / "singularity_cache"
    cache_path.mkdir()
    with monkeypatch.context() as m:
        m.setenv("PATH", new_paths)
        m.setenv("CWL_SINGULARITY_CACHE", str(cache_path))
        js_engine = sandboxjs.get_js_engine()
        js_engine.localdata = threading.local()  # type: ignore[attr-defined]
        js_engine.have_node_slim = False  # type: ignore[attr-defined]
        assert (
            expression.do_eval(
                "$(42*23)",
                {},
                [{"class": "InlineJavascriptRequirement"}],
                None,
                None,
                {},
                cwlVersion="v1.0",
                container_engine="singularity",
            )
            == 42 * 23
        )
        assert (cache_path / "node_alpine.sif").exists()


def test_caches_js_processes(mocker: Any) -> None:
    sandboxjs.exec_js_process("7", context="{}")

    mocked_new_js_proc = mocker.patch("cwl_utils.sandboxjs.new_js_proc")
    sandboxjs.exec_js_process("7", context="{}")

    mocked_new_js_proc.assert_not_called()
