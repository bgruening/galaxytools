"""
An interface the clients should implement.

This class is primarily a helper for the library and user code
should not use it directly.
"""

import time
from typing import (
    Any,
    Literal,
    Optional,
    overload,
    TYPE_CHECKING,
)

import requests

import bioblend

# The following import must be preserved for compatibility because
# ConnectionError class was originally defined here
from bioblend import ConnectionError

if TYPE_CHECKING:
    from bioblend.galaxyclient import GalaxyClient


class Client:
    # The `module` attribute needs to be defined in subclasses
    module: str
    gi: "GalaxyClient"

    @classmethod
    def max_get_retries(cls) -> None:
        raise AttributeError("Deprecated method, please use gi's `max_get_attempts` property")

    @classmethod
    def set_max_get_retries(cls, value: int) -> None:
        raise AttributeError("Deprecated method, please use gi's `max_get_attempts` property")

    @classmethod
    def get_retry_delay(cls) -> None:
        raise AttributeError("Deprecated method, please use gi's `get_retry_delay` property")

    @classmethod
    def set_get_retry_delay(cls, value: float) -> None:
        raise AttributeError("Deprecated method, please use gi's `get_retry_delay` property")

    def __init__(self, galaxy_instance: "GalaxyClient") -> None:
        """
        A generic Client interface defining the common fields.

        All clients *must* define the following field (which will be
        used as part of the URL composition (e.g.,
        ``http://<galaxy_instance>/api/libraries``): ``self.module =
        'workflows' | 'libraries' | 'histories' | ...``
        """
        self.gi = galaxy_instance

    def _make_url(self, module_id: Optional[str] = None, deleted: bool = False, contents: bool = False) -> str:
        """
        Compose a URL based on the provided arguments.

        :type module_id: str
        :param module_id: The encoded ID for a specific module (eg, library ID)

        :type deleted: bool
        :param deleted: If ``True``, include ``deleted`` in the URL, after the module
                        name (eg, ``<base_url>/api/libraries/deleted``)

        :type contents: bool
        :param contents: If ``True``, include 'contents' in the URL, after the module ID:
                         ``<base_url>/api/libraries/<encoded_library_id>/contents``
        """
        c_url = "/".join((self.gi.url, self.module))
        if deleted:
            c_url = c_url + "/deleted"
        if module_id:
            c_url = "/".join((c_url, module_id))
            if contents:
                c_url = c_url + "/contents"
        return c_url

    @overload
    def _get(
        self,
        id: Optional[str] = None,
        deleted: bool = False,
        contents: bool = False,
        url: Optional[str] = None,
        params: Optional[dict] = None,
        *,
        json: Literal[False],
        stream: bool = False,
    ) -> requests.Response: ...

    @overload
    def _get(
        self,
        id: Optional[str] = None,
        deleted: bool = False,
        contents: bool = False,
        url: Optional[str] = None,
        params: Optional[dict] = None,
        json: bool = True,
        stream: bool = False,
    ) -> Any: ...

    def _get(
        self,
        id: Optional[str] = None,
        deleted: bool = False,
        contents: bool = False,
        url: Optional[str] = None,
        params: Optional[dict] = None,
        json: bool = True,
        stream: bool = False,
    ) -> Any:
        """
        Do a GET request, composing the URL from ``id``, ``deleted`` and
        ``contents``.  Alternatively, an explicit ``url`` can be provided.
        If ``json`` is set to ``True``, return a decoded JSON object
        (and treat an empty or undecodable response as an error).

        The request will optionally be retried as configured by gi's
        ``max_get_attempts`` and ``get_retry_delay``: this offers some
        resilience in the presence of temporary failures.

        :return: The decoded response if ``json`` is set to ``True``, otherwise
          the response object
        """
        if url is None:
            url = self._make_url(module_id=id, deleted=deleted, contents=contents)
        attempts_left = self.gi.max_get_attempts
        retry_delay = self.gi.get_retry_delay
        bioblend.log.debug("GET - attempts left: %s; retry delay: %s", attempts_left, retry_delay)
        msg = ""
        while attempts_left > 0:
            attempts_left -= 1
            try:
                r = self.gi.make_get_request(url, params=params, stream=stream)
            except requests.exceptions.ConnectionError as e:
                msg = str(e)
                r = requests.Response()  # empty Response object used when raising ConnectionError
            else:
                if r.status_code == 200:
                    if not json:
                        return r
                    elif not r.content:
                        msg = "GET: empty response"
                    else:
                        try:
                            return r.json()
                        except ValueError:
                            msg = f"GET: invalid JSON : {r.content!r}"
                else:
                    msg = f"GET: error {r.status_code}: {r.content!r}"
            msg = f"{msg}, {attempts_left} attempts left"
            if attempts_left <= 0:
                bioblend.log.error(msg)
                raise ConnectionError(
                    msg,
                    body=r.text,
                    status_code=r.status_code,
                )
            else:
                bioblend.log.warning(msg)
                time.sleep(retry_delay)

    def _post(
        self,
        payload: Optional[dict] = None,
        id: Optional[str] = None,
        deleted: bool = False,
        contents: bool = False,
        url: Optional[str] = None,
        files_attached: bool = False,
    ) -> Any:
        """
        Do a generic POST request, composing the url from the contents of the
        arguments. Alternatively, an explicit ``url`` can be provided to use
        for the request.
        The payload dict may contain file handles (in which case the
        ``files_attached`` flag must be set to true).

        If ``files_attached`` is set to ``False``, the request body will be
        JSON-encoded; otherwise, it will be encoded as multipart/form-data.

        :type payload: dict
        :param payload: additional parameters to send in the body of the request

        :return: The decoded response.
        """
        if not url:
            url = self._make_url(module_id=id, deleted=deleted, contents=contents)
        return self.gi.make_post_request(url, payload=payload, files_attached=files_attached)

    def _put(
        self,
        payload: Optional[dict] = None,
        id: Optional[str] = None,
        url: Optional[str] = None,
        params: Optional[dict] = None,
    ) -> Any:
        """
        Do a generic PUT request, composing the url from the contents of the
        arguments. Alternatively, an explicit ``url`` can be provided to use
        for the request.

        :type payload: dict
        :param payload: additional parameters to send in the body of the request

        :return: The decoded response.
        """
        if not url:
            url = self._make_url(module_id=id)
        return self.gi.make_put_request(url, payload=payload, params=params)

    def _patch(
        self,
        payload: Optional[dict] = None,
        id: Optional[str] = None,
        url: Optional[str] = None,
        params: Optional[dict] = None,
    ) -> Any:
        """
        Do a generic PATCH request, composing the url from the contents of the
        arguments. Alternatively, an explicit ``url`` can be provided to use
        for the request.

        :type payload: dict
        :param payload: additional parameters to send in the body of the request

        :return: The decoded response.
        """
        if not url:
            url = self._make_url(module_id=id)
        return self.gi.make_patch_request(url, payload=payload, params=params)

    def _delete(
        self,
        payload: Optional[dict] = None,
        id: Optional[str] = None,
        deleted: bool = False,
        contents: bool = False,
        url: Optional[str] = None,
        params: Optional[dict] = None,
    ) -> Any:
        """
        Do a generic DELETE request, composing the url from the contents of the
        arguments. Alternatively, an explicit ``url`` can be provided to use
        for the request.

        :type payload: dict
        :param payload: additional parameters to send in the body of the request

        :return: The decoded response or None.
        """
        if not url:
            url = self._make_url(module_id=id, deleted=deleted, contents=contents)
        r = self.gi.make_delete_request(url, payload=payload, params=params)
        if 200 <= r.status_code < 202:
            return r.json()
        elif 202 <= r.status_code < 300:
            return None
        # @see self.body for HTTP response body
        raise ConnectionError(
            f"Unexpected HTTP status code: {r.status_code}",
            body=r.text,
            status_code=r.status_code,
        )
