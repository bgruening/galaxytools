"""Resource fetching."""

import logging
import os
import re
import sys
import urllib.parse
import urllib.request
from abc import ABC, abstractmethod
from typing import Final, Optional

import requests
from mypy_extensions import mypyc_attr

from .exceptions import ValidationException
from .utils import CacheType

_re_drive: Final = re.compile(r"/([a-zA-Z]):")
_logger: Final = logging.getLogger("salad")


@mypyc_attr(allow_interpreted_subclasses=True)
class Fetcher(ABC):
    """Fetch resources from URIs."""

    @abstractmethod
    def fetch_text(self, url: str, content_types: Optional[list[str]] = None) -> str:
        """Retrieve the given resource as a string."""

    @abstractmethod
    def check_exists(self, url: str) -> bool:
        """Check if the given resource exists."""

    @abstractmethod
    def urljoin(self, base_url: str, url: str) -> str:
        """Construct a full (“absolute”) URL by combining a “base URL” with another URL."""

    schemes = ["file", "http", "https", "mailto"]

    def supported_schemes(self) -> list[str]:
        """Return the list of supported URI schemes."""
        return self.schemes


@mypyc_attr(allow_interpreted_subclasses=True)
class MemoryCachingFetcher(Fetcher):
    """Fetcher that caches resources in memory after retrieval."""

    def __init__(self, cache: CacheType) -> None:
        """Create a MemoryCachingFetcher object."""
        self.cache = cache


@mypyc_attr(allow_interpreted_subclasses=True)
class DefaultFetcher(MemoryCachingFetcher):
    """The default Fetcher implementation."""

    def __init__(
        self,
        cache: CacheType,
        session: Optional[requests.sessions.Session],
    ) -> None:
        """Create a DefaultFetcher object."""
        super().__init__(cache)
        self.session: Final = session

    def fetch_text(self, url: str, content_types: Optional[list[str]] = None) -> str:
        """Retrieve the given resource as a string."""
        result: Final = self.cache.get(url, None)
        if isinstance(result, str):
            return result

        split: Final = urllib.parse.urlsplit(url)
        scheme: Final = split.scheme
        path = split.path

        if scheme in ["http", "https"] and self.session is not None:
            try:
                headers = {}
                if content_types:
                    headers["Accept"] = ", ".join(content_types) + ", */*;q=0.8"
                resp: Final = self.session.get(url, headers=headers)
                resp.raise_for_status()
            except Exception as e:
                raise ValidationException(f"Error fetching {url}: {e}") from e
            if content_types and "content-type" in resp.headers:
                received_content_types = set(
                    resp.headers["content-type"].split(";")[:1][0].split(",")
                )
                if set(content_types).isdisjoint(received_content_types):
                    _logger.warning(
                        "While fetching %s, got content-type of %r. Expected one of %s.",
                        url,
                        resp.headers["content-type"].split(";")[:1][0],
                        content_types,
                    )
            return resp.text
        if scheme == "file":
            try:
                # On Windows, url.path will be /drive:/path ; on Unix systems,
                # /path. As we want drive:/path instead of /drive:/path on Windows,
                # remove the leading /.
                if os.path.isabs(
                    path[1:]
                ):  # checking if pathis valid after removing front / or not
                    path = path[1:]
                with open(urllib.request.url2pathname(str(path)), encoding="utf-8") as fp:
                    return str(fp.read())

            except OSError as err:
                if err.filename == path:
                    raise ValidationException(str(err)) from err
                raise ValidationException(f"Error reading {url}: {err}") from err
        raise ValidationException(f"Unsupported scheme in url: {url}")

    def check_exists(self, url: str) -> bool:
        if url in self.cache:
            return True

        split: Final = urllib.parse.urlsplit(url)
        scheme: Final = split.scheme
        path: Final = split.path

        if scheme in ["http", "https"]:
            if self.session is None:
                raise ValidationException(f"Can't check {scheme} URL, session is None")
            try:
                resp: Final = self.session.head(url, allow_redirects=True)
                resp.raise_for_status()
            except Exception:
                return False
            self.cache[url] = True
            return True
        if scheme == "file":
            return os.path.exists(urllib.request.url2pathname(str(path)))
        if scheme == "mailto":
            return True
        raise ValidationException(f"Unsupported scheme {scheme!r} in url: {url}")

    def urljoin(self, base_url: str, url: str) -> str:
        if url.startswith("_:"):
            return url

        basesplit = urllib.parse.urlsplit(base_url)
        split = urllib.parse.urlsplit(url)
        if basesplit.scheme and basesplit.scheme != "file" and split.scheme == "file":
            raise ValidationException(
                f"Not resolving potential remote exploit {url} from base {base_url}"
            )

        if sys.platform == "win32":
            if base_url == url:
                return url
            basesplit = urllib.parse.urlsplit(base_url)
            # note that below might split
            # "C:" with "C" as URI scheme
            split = urllib.parse.urlsplit(url)

            has_drive = split.scheme and len(split.scheme) == 1

            if basesplit.scheme == "file":
                # Special handling of relative file references on Windows
                # as urllib seems to not be quite up to the job

                # netloc MIGHT appear in equivalents of UNC Strings
                # \\server1.example.com\path as
                # file:///server1.example.com/path
                # https://tools.ietf.org/html/rfc8089#appendix-E.3.2
                # (TODO: test this)
                netloc: Final = split.netloc or basesplit.netloc

                # Check if url is a local path like "C:/Users/fred"
                # or actually an absolute URI like http://example.com/fred
                if has_drive:
                    # Assume split.scheme is actually a drive, e.g. "C:"
                    # so we'll recombine into a path
                    path_with_drive = urllib.parse.urlunsplit(
                        (split.scheme, "", split.path, "", "")
                    )
                    # Compose new file:/// URI with path_with_drive
                    # .. carrying over any #fragment (?query just in case..)
                    return urllib.parse.urlunsplit(
                        ("file", netloc, path_with_drive, split.query, split.fragment)
                    )
                if not split.scheme and not netloc and split.path and split.path.startswith("/"):
                    # Relative - but does it have a drive?
                    base_drive: Final = _re_drive.match(basesplit.path)
                    drive: Final = _re_drive.match(split.path)
                    if base_drive and not drive:
                        # Keep drive letter from base_url
                        # https://tools.ietf.org/html/rfc8089#appendix-E.2.1
                        # e.g. urljoin("file:///D:/bar/a.txt", "/foo/b.txt")
                        #          == file:///D:/foo/b.txt
                        path_with_drive = f"/{base_drive.group(1)}:{split.path}"
                        return urllib.parse.urlunsplit(
                            (
                                "file",
                                netloc,
                                path_with_drive,
                                split.query,
                                split.fragment,
                            )
                        )

                # else: fall-through to resolve as relative URI
            elif has_drive:
                # Base is http://something but url is C:/something - which urllib
                # would wrongly resolve as an absolute path that could later be used
                # to access local files
                raise ValidationException(
                    f"Not resolving potential remote exploit {url} from base {base_url}"
                )

        return urllib.parse.urljoin(base_url, url)
