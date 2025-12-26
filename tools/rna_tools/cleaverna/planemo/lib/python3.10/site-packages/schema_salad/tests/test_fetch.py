import os
from typing import Optional
from urllib.parse import urljoin, urlsplit

import pytest
import requests

from schema_salad.fetcher import Fetcher
from schema_salad.ref_resolver import Loader, file_uri
from schema_salad.utils import CacheType


class testFetcher(Fetcher):
    def __init__(
        self,
        cache: CacheType,
        session: Optional[requests.sessions.Session],
    ) -> None:
        pass

    def fetch_text(self, url: str, content_types: Optional[list[str]] = None) -> str:
        if url == "keep:abc+123/foo.txt":
            return "hello: keepfoo"
        if url.endswith("foo.txt"):
            return "hello: foo"
        else:
            raise RuntimeError("Not foo.txt")

    def check_exists(self, url: str) -> bool:
        if url.endswith("foo.txt"):
            return True
        else:
            return False

    def urljoin(self, base: str, url: str) -> str:
        urlsp = urlsplit(url)
        if urlsp.scheme:
            return url
        basesp = urlsplit(base)

        if basesp.scheme == "keep":
            return base + "/" + url
        return urljoin(base, url)


class CWLTestFetcher(Fetcher):
    def __init__(
        self,
        cache: CacheType,
        session: Optional[requests.sessions.Session],
    ) -> None:
        pass

    def fetch_text(self, url: str, content_types: Optional[list[str]] = None) -> str:
        if url == "baz:bar/foo.cwl":
            return """
cwlVersion: v1.0
class: CommandLineTool
baseCommand: echo
inputs: []
outputs: []
"""
        raise RuntimeError(f"Not foo.cwl, was {url}")

    def check_exists(self, url: str) -> bool:
        return url == "baz:bar/foo.cwl"

    def urljoin(self, base: str, url: str) -> str:
        urlsp = urlsplit(url)
        if urlsp.scheme:
            return url
        basesp = urlsplit(base)

        if basesp.scheme == "keep":
            return base + "/" + url
        return urljoin(base, url)


def test_fetcher() -> None:
    loader = Loader({}, fetcher_constructor=testFetcher)
    assert {"hello": "foo"} == loader.resolve_ref("foo.txt")[0]
    assert {"hello": "keepfoo"} == loader.resolve_ref("foo.txt", base_url="keep:abc+123")[0]
    assert loader.check_exists("foo.txt")

    with pytest.raises(RuntimeError):
        loader.resolve_ref("bar.txt")
    assert not loader.check_exists("bar.txt")


def test_cache() -> None:
    loader = Loader({})
    foo = os.path.join(os.getcwd(), "foo.txt")
    foo = file_uri(foo)
    loader.cache.update({foo: "hello: foo"})
    print(loader.cache)
    assert {"hello": "foo"} == loader.resolve_ref("foo.txt")[0]
    assert loader.check_exists(foo)
