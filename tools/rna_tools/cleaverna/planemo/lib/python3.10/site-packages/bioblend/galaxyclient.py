"""
Helper class for Galaxy and ToolShed Instance object

This class is primarily a helper for the library and user code
should not use it directly.
A base representation of an instance
"""

import base64
import contextlib
import json
import logging
from typing import (
    Any,
    Optional,
    Union,
)

import requests
import tusclient.client
import tusclient.exceptions
import tusclient.storage.filestorage
import tusclient.uploader
from requests_toolbelt import MultipartEncoder

from bioblend import ConnectionError
from bioblend.util import FileStream

log = logging.getLogger(__name__)

UPLOAD_CHUNK_SIZE = 10**7


class GalaxyClient:
    def __init__(
        self,
        url: str,
        key: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None,
        *,
        token: Optional[str] = None,
        verify: bool = True,
        timeout: Optional[float] = None,
        user_agent: Optional[str] = None,
    ) -> None:
        """
        :param verify: Whether to verify the server's TLS certificate
        :type verify: bool
        :param timeout: Timeout for requests operations, set to None for no timeout (the default).
        :type timeout: float
        """
        self.verify = verify
        self.timeout = timeout
        # Make sure the URL scheme is defined (otherwise requests will not work)
        if not url.lower().startswith("http"):
            found_scheme = None
            # Try to guess the scheme, starting from the more secure
            for scheme in ("https://", "http://"):
                log.warning("Missing scheme in url, trying with %s", scheme)
                with contextlib.suppress(requests.RequestException):
                    r = requests.get(
                        scheme + url,
                        timeout=self.timeout,
                        verify=self.verify,
                    )
                    r.raise_for_status()
                    found_scheme = scheme
                    break
            else:
                raise ValueError(f"Missing scheme in url {url}")
            url = found_scheme + url
        self.base_url = url.rstrip("/")
        # All of Galaxy's and ToolShed's API's are rooted at <url>/api so make that the url
        self.url = f"{self.base_url}/api"
        # If key has been supplied, use it; otherwise just set email and
        # password and grab user's key before first request.
        if key:
            self._key: Optional[str] = key
        elif token:
            self.token: Optional[str] = token
        else:
            self._key = None
            self.email = email
            self.password = password
        self.json_headers: dict[str, Union[str, bytes, None]] = {"Content-Type": "application/json"}
        if user_agent:
            self.json_headers["User-Agent"] = user_agent
        # json_headers needs to be set before key can be defined, otherwise authentication with email/password causes an error
        if token:
            self.json_headers["Authorization"] = f"Bearer {token}"
        else:
            self.json_headers["x-api-key"] = self.key
        # Number of attempts before giving up on a GET request.
        self._max_get_attempts = 1
        # Delay in seconds between subsequent retries.
        self._get_retry_delay = 10.0

    @property
    def max_get_attempts(self) -> int:
        """
        The maximum number of attempts for a GET request. Default: 1
        """
        return self._max_get_attempts

    @max_get_attempts.setter
    def max_get_attempts(self, value: int) -> None:
        """
        Set the maximum number of attempts for GET requests. A value greater
        than one causes failed GET requests to be retried `value` - 1 times.
        """
        if value < 1:
            raise ValueError(f"Number of attempts must be >= 1 (got: {value})")
        self._max_get_attempts = value

    @property
    def get_retry_delay(self) -> float:
        """
        The delay (in seconds) to wait before retrying a failed GET request.
        Default: 10.0
        """
        return self._get_retry_delay

    @get_retry_delay.setter
    def get_retry_delay(self, value: float) -> None:
        """
        Set the delay (in seconds) to wait before retrying a failed GET
        request.
        """
        if value < 0:
            raise ValueError(f"Retry delay must be >= 0 (got: {value})")
        self._get_retry_delay = value

    def make_get_request(self, url: str, **kwargs: Any) -> requests.Response:
        """
        Make a GET request using the provided ``url``.

        Keyword arguments are the same as in requests.request.

        If ``verify`` is not provided, ``self.verify`` will be used.

        :rtype: requests.Response
        :return: the response object.
        """
        headers = self.json_headers
        kwargs.setdefault("timeout", self.timeout)
        kwargs.setdefault("verify", self.verify)
        r = requests.get(url, headers=headers, **kwargs)
        return r

    def make_post_request(
        self, url: str, payload: Optional[dict] = None, params: Optional[dict] = None, files_attached: bool = False
    ) -> Any:
        """
        Make a POST request using the provided ``url`` and ``payload``.
        The ``payload`` must be a dict that contains the request values.
        The payload dict may contain file handles (in which case the files_attached
        flag must be set to true).

        :return: The decoded response.
        """

        def my_dumps(d: dict) -> dict:
            """
            Apply ``json.dumps()`` to the values of the dict ``d`` if they are
            not of type ``FileStream``.
            """
            for k, v in d.items():
                if not isinstance(v, (FileStream, str, bytes)):
                    d[k] = json.dumps(v)
            return d

        # Compute data, headers, params arguments for request.post,
        # leveraging the requests-toolbelt library if any files have
        # been attached.
        if files_attached:
            payload_copy = payload.copy() if payload is not None else {}
            if params:
                payload_copy.update(params)
            data = MultipartEncoder(fields=my_dumps(payload_copy))
            headers = self.json_headers.copy()
            headers["Content-Type"] = data.content_type
            post_params = None
        else:
            data = json.dumps(payload) if payload is not None else None
            headers = self.json_headers
            post_params = params

        r = requests.post(
            url,
            params=post_params,
            data=data,
            headers=headers,
            timeout=self.timeout,
            allow_redirects=False,
            verify=self.verify,
        )
        if r.status_code == 200:
            try:
                return r.json()
            except Exception as e:
                raise ConnectionError(
                    f"Request was successful, but cannot decode the response content: {e}",
                    body=r.content,
                    status_code=r.status_code,
                )
        # @see self.body for HTTP response body
        raise ConnectionError(
            f"Unexpected HTTP status code: {r.status_code}",
            body=r.text,
            status_code=r.status_code,
        )

    def make_delete_request(
        self, url: str, payload: Optional[dict] = None, params: Optional[dict] = None
    ) -> requests.Response:
        """
        Make a DELETE request using the provided ``url`` and the optional
        arguments.

        :type payload: dict
        :param payload: a JSON-serializable dictionary

        :rtype: requests.Response
        :return: the response object.
        """
        data = json.dumps(payload) if payload is not None else None
        headers = self.json_headers
        r = requests.delete(
            url,
            params=params,
            data=data,
            headers=headers,
            timeout=self.timeout,
            allow_redirects=False,
            verify=self.verify,
        )
        return r

    def make_put_request(self, url: str, payload: Optional[dict] = None, params: Optional[dict] = None) -> Any:
        """
        Make a PUT request using the provided ``url`` with required payload.

        :type payload: dict
        :param payload: a JSON-serializable dictionary

        :return: The decoded response.
        """
        data = json.dumps(payload) if payload is not None else None
        headers = self.json_headers
        r = requests.put(
            url,
            params=params,
            data=data,
            headers=headers,
            timeout=self.timeout,
            allow_redirects=False,
            verify=self.verify,
        )
        if r.status_code == 200:
            try:
                return r.json()
            except Exception as e:
                raise ConnectionError(
                    f"Request was successful, but cannot decode the response content: {e}",
                    body=r.content,
                    status_code=r.status_code,
                )
        # @see self.body for HTTP response body
        raise ConnectionError(
            f"Unexpected HTTP status code: {r.status_code}",
            body=r.text,
            status_code=r.status_code,
        )

    def make_patch_request(self, url: str, payload: Optional[dict] = None, params: Optional[dict] = None) -> Any:
        """
        Make a PATCH request using the provided ``url`` with required payload.

        :type payload: dict
        :param payload: a JSON-serializable dictionary

        :return: The decoded response.
        """
        data = json.dumps(payload) if payload is not None else None
        headers = self.json_headers
        r = requests.patch(
            url,
            params=params,
            data=data,
            headers=headers,
            timeout=self.timeout,
            allow_redirects=False,
            verify=self.verify,
        )
        if r.status_code == 200:
            try:
                return r.json()
            except Exception as e:
                raise ConnectionError(
                    f"Request was successful, but cannot decode the response content: {e}",
                    body=r.content,
                    status_code=r.status_code,
                )
        # @see self.body for HTTP response body
        raise ConnectionError(
            f"Unexpected HTTP status code: {r.status_code}",
            body=r.text,
            status_code=r.status_code,
        )

    def get_tus_uploader(
        self,
        path: str,
        url: str = "/upload/resumable_upload",
        storage: Optional[str] = None,
        metadata: Optional[dict] = None,
        chunk_size: Optional[int] = UPLOAD_CHUNK_SIZE,
    ) -> tusclient.uploader.Uploader:
        """
        Return the tus client uploader object for uploading to the Galaxy tus endpoint

        :type path: str
        :param path: path of the file to upload

        :type url: str
        :param url: URL (relative to base URL) of the upload endpoint

        :type storage: str
        :param storage: Local path to store URLs resuming uploads

        :type metadata: dict
        :param metadata: Metadata to send with upload request

        :type chunk_size: int
        :param chunk_size: Number of bytes to send in each chunk

        :rtype: tusclient.uploader.Uploader
        :return: tus uploader object
        """
        headers = {"x-api-key": self.key}
        client = tusclient.client.TusClient(self.url + url, headers=headers)
        url_storage = tusclient.storage.filestorage.FileStorage(storage) if storage else None
        try:
            return client.uploader(
                file_path=path,
                chunk_size=chunk_size,
                metadata=metadata,
                store_url=url_storage is not None,
                url_storage=url_storage,
            )
        except tusclient.exceptions.TusCommunicationError as exc:
            raise ConnectionError(
                f"Unexpected HTTP status code: {exc.status_code}",
                body=str(exc),
                status_code=exc.status_code,
            )

    @property
    def key(self) -> Optional[str]:
        if not self._key and self.email is not None and self.password is not None:
            unencoded_credentials = f"{self.email}:{self.password}"
            authorization = base64.b64encode(unencoded_credentials.encode())
            headers = self.json_headers.copy()
            headers["Authorization"] = authorization
            auth_url = f"{self.url}/authenticate/baseauth"
            # Use lower level method instead of make_get_request() because we
            # need the additional Authorization header.
            r = requests.get(
                auth_url,
                headers=headers,
                timeout=self.timeout,
                verify=self.verify,
            )
            if r.status_code != 200:
                raise Exception("Failed to authenticate user.")
            response = r.json()
            if isinstance(response, str):
                # bug in Tool Shed
                response = json.loads(response)
            self._key = response["api_key"]
        return self._key


def _tus_uploader_session_id(self: tusclient.uploader.Uploader) -> str:
    assert self.url
    return self.url.rsplit("/", 1)[1]


# monkeypatch a session_id property on to uploader
tusclient.uploader.Uploader.session_id = property(_tus_uploader_session_id)
