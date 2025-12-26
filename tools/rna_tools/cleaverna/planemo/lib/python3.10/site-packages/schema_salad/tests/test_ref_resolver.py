"""Test the ref_resolver module."""

import os
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any, Union

import pytest
from _pytest.fixtures import FixtureRequest
from requests import Session
from ruamel.yaml.comments import CommentedMap, CommentedSeq

import schema_salad.main
from schema_salad.exceptions import ValidationException
from schema_salad.fetcher import DefaultFetcher
from schema_salad.ref_resolver import Loader
from schema_salad.tests.util import get_data, get_data_uri, get_path


def is_fs_case_sensitive(
    path: Path,
) -> bool:  # https://stackoverflow.com/a/36612604/1585509
    with tempfile.NamedTemporaryFile(prefix="TmP", dir=path) as tmp_file:
        return not os.path.exists(tmp_file.name.lower())


@pytest.fixture
def tmp_dir_fixture(request: FixtureRequest) -> str:
    d = tempfile.mkdtemp()

    @request.addfinalizer
    def teardown() -> None:
        shutil.rmtree(d)

    return d


def test_Loader_initialisation_for_HOME_env_var(tmp_dir_fixture: str) -> None:
    # Ensure HOME is set.
    os.environ["HOME"] = tmp_dir_fixture

    loader = Loader(ctx={})
    assert isinstance(loader.session, Session)


def test_Loader_initialisation_for_TMP_env_var(tmp_dir_fixture: str) -> None:
    # Ensure HOME is missing.
    if "HOME" in os.environ:
        del os.environ["HOME"]
    # Ensure TMP is present.
    os.environ["TMP"] = tmp_dir_fixture

    loader = Loader(ctx={})
    assert isinstance(loader.session, Session)


def test_Loader_initialisation_with_neither_TMP_HOME_set(tmp_dir_fixture: str) -> None:
    # Ensure HOME is missing.
    if "HOME" in os.environ:
        del os.environ["HOME"]
    if "TMP" in os.environ:
        del os.environ["TMP"]

    loader = Loader(ctx={})
    assert isinstance(loader.session, Session)


def test_Loader_initialisation_disable_doc_cache(tmp_dir_fixture: str) -> None:
    loader = Loader(ctx={}, doc_cache=False)
    assert isinstance(loader.session, Session)


@pytest.mark.skipif(sys.platform != "win32", reason="Only for win32")
def test_DefaultFetcher_urljoin_win32(tmp_dir_fixture: str) -> None:
    # Ensure HOME is set.
    os.environ["HOME"] = tmp_dir_fixture

    fetcher = DefaultFetcher({}, None)
    # Relative path, same folder
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "soup.cwl")
    assert url == "file:///C:/Users/fred/soup.cwl"
    # Relative path, sub folder
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "foo/soup.cwl")
    assert url == "file:///C:/Users/fred/foo/soup.cwl"
    # relative climb-up path
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "../alice/soup.cwl")
    assert url == "file:///C:/Users/alice/soup.cwl"

    # Path with drive: should not be treated as relative to directory
    # Note: \ would already have been converted to / by resolve_ref()
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "c:/bar/soup.cwl")
    assert url == "file:///c:/bar/soup.cwl"
    # /C:/  (regular URI absolute path)
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "/c:/bar/soup.cwl")
    assert url == "file:///c:/bar/soup.cwl"
    # Relative, change drive
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "D:/baz/soup.cwl")
    assert url == "file:///d:/baz/soup.cwl"
    # Relative from root of base's D: drive
    url = fetcher.urljoin("file:///d:/baz/soup.cwl", "/foo/soup.cwl")
    assert url == "file:///d:/foo/soup.cwl"

    # resolving absolute non-drive URIs still works
    url = fetcher.urljoin("file:///C:/Users/fred/foo.cwl", "http://example.com/bar/soup.cwl")
    assert url == "http://example.com/bar/soup.cwl"
    # and of course relative paths from http://
    url = fetcher.urljoin("http://example.com/fred/foo.cwl", "soup.cwl")
    assert url == "http://example.com/fred/soup.cwl"

    # Stay on http:// and same host
    url = fetcher.urljoin("http://example.com/fred/foo.cwl", "/bar/soup.cwl")
    assert url == "http://example.com/bar/soup.cwl"

    # Security concern - can't resolve file: from http:
    with pytest.raises(ValidationException):
        url = fetcher.urljoin("http://example.com/fred/foo.cwl", "file:///c:/bar/soup.cwl")
    # Drive-relative -- should NOT return "absolute" URI c:/bar/soup.cwl"
    # as that is a potential remote exploit
    with pytest.raises(ValidationException):
        url = fetcher.urljoin("http://example.com/fred/foo.cwl", "c:/bar/soup.cwl")


def test_DefaultFetcher_urljoin_linux(tmp_dir_fixture: str) -> None:
    # Ensure HOME is set.
    os.environ["HOME"] = tmp_dir_fixture

    actual_platform = sys.platform
    try:
        # Pretend it's Linux (e.g. not win32)
        sys.platform = "linux2"
        fetcher = DefaultFetcher({}, None)
        url = fetcher.urljoin("file:///home/fred/foo.cwl", "soup.cwl")
        assert url == "file:///home/fred/soup.cwl"

        url = fetcher.urljoin("file:///home/fred/foo.cwl", "../alice/soup.cwl")
        assert url == "file:///home/alice/soup.cwl"
        # relative from root
        url = fetcher.urljoin("file:///home/fred/foo.cwl", "/baz/soup.cwl")
        assert url == "file:///baz/soup.cwl"

        url = fetcher.urljoin("file:///home/fred/foo.cwl", "http://example.com/bar/soup.cwl")
        assert url == "http://example.com/bar/soup.cwl"

        url = fetcher.urljoin("http://example.com/fred/foo.cwl", "soup.cwl")
        assert url == "http://example.com/fred/soup.cwl"

        # Root-relative -- here relative to http host, not file:///
        url = fetcher.urljoin("http://example.com/fred/foo.cwl", "/bar/soup.cwl")
        assert url == "http://example.com/bar/soup.cwl"

        # Security concern - can't resolve file: from http:
        with pytest.raises(ValidationException):
            url = fetcher.urljoin("http://example.com/fred/foo.cwl", "file:///bar/soup.cwl")

        # But this one is not "dangerous" on Linux
        fetcher.urljoin("http://example.com/fred/foo.cwl", "c:/bar/soup.cwl")

    finally:
        sys.platform = actual_platform


def test_import_list() -> None:
    import schema_salad.ref_resolver
    from schema_salad.sourceline import cmap

    basedir = schema_salad.ref_resolver.file_uri(os.path.dirname(__file__) + "/")
    loader = schema_salad.ref_resolver.Loader({})
    ra, _ = loader.resolve_all(cmap({"foo": {"$import": "list.json"}}), basedir)

    assert {"foo": ["bar", "baz"]} == ra

    importfile = "import:%slist.json" % basedir

    assert importfile in loader.idx


def test_fetch_inject_id() -> None:
    path = get_path("tests/inject-id1.yml")
    if is_fs_case_sensitive(path.parent):

        def lower(item: str) -> str:
            return item

    else:

        def lower(item: str) -> str:
            return item.lower()

    l1 = Loader({"id": "@id"})
    furi1 = path.as_uri()
    r1, _ = l1.resolve_ref(furi1)
    assert {"id": furi1 + "#foo", "bar": "baz"} == r1
    assert [lower(furi1), lower(furi1 + "#foo")] == sorted(list(lower(k) for k in l1.idx.keys()))

    l2 = Loader({"id": "@id"})
    path2_uri = get_data_uri("tests/inject-id2.yml")
    r2, _ = l2.resolve_ref(path2_uri)
    assert {"id": path2_uri, "bar": "baz"} == r2
    assert [lower(path2_uri)] == sorted(list(lower(k) for k in l2.idx.keys()))

    l3 = Loader({"id": "@id"})
    path3_uri = get_data_uri("tests/inject-id3.yml")
    r3, _ = l3.resolve_ref(path3_uri)
    assert {"id": "http://example.com", "bar": "baz"} == r3
    assert [lower(path3_uri), "http://example.com"] == sorted(list(lower(k) for k in l3.idx.keys()))


def test_attachments() -> None:
    path = get_path("tests/multidoc.yml")

    furi = path.as_uri()

    l1 = Loader({})
    r1, _ = l1.resolve_ref(furi)
    with path.open() as f1:
        content = f1.read()
        assert {"foo": "bar", "baz": content, "quux": content} == r1

    def aa1(item: Union[CommentedMap, CommentedSeq]) -> bool:
        return bool(item["foo"] == "bar")

    l2 = Loader({}, allow_attachments=aa1)
    r2, _ = l2.resolve_ref(furi)
    assert {
        "foo": "bar",
        "baz": "This is the {first attachment}.\n",
        "quux": "This is the [second attachment].",
    } == r2

    def aa2(item: Union[CommentedMap, CommentedSeq]) -> bool:
        return bool(item["foo"] == "baz")

    l3 = Loader({}, allow_attachments=aa2)
    r3, _ = l3.resolve_ref(furi)
    with path.open() as f2:
        content = f2.read()
        assert {"foo": "bar", "baz": content, "quux": content} == r3


def test_check_exists_follows_redirects() -> None:
    fetcher = DefaultFetcher({}, Session())
    # This tests that the redirect from http to https doesn't cause a
    # false positive, this URL doesn't exist and it should return
    # False because if it follows the redirect it'll get a 404.  Tests
    # for previous bug where allow_redirects wasn't set and as a
    # result it would return True even though the result was 3xx.
    assert not fetcher.check_exists("http://commonwl.org/does/not/exist")


def test_resolve_missing_step_id(caplog: Any) -> None:
    """From issue #cwltool/issues/1635. A Workflow with a Step without
    the name attribute must raise a ValidationException that contains
    the SourceLine data."""
    schema_path = get_data("tests/test_schema/CommonWorkflowLanguage.yml")
    document_path = get_data("tests/missing_step_name.cwl")
    assert 1 == schema_salad.main.main(argsl=["--print-rdf", schema_path, document_path])
    assert "missing_step_name.cwl:13:1" in "\n".join(caplog.messages)
