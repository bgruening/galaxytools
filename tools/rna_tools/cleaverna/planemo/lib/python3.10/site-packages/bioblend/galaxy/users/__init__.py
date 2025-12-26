"""
Contains possible interaction dealing with Galaxy users.

Most of these methods must be executed by a registered Galaxy admin user.
"""

from typing import (
    Any,
    Optional,
    TYPE_CHECKING,
)

from bioblend import ConnectionError
from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class UserClient(Client):
    module = "users"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_users(
        self,
        deleted: bool = False,
        f_email: Optional[str] = None,
        f_name: Optional[str] = None,
        f_any: Optional[str] = None,
    ) -> list[dict[str, Any]]:
        """
        Get a list of all registered users. If ``deleted`` is set to ``True``,
        get a list of deleted users.

        :type deleted: bool
        :param deleted: Whether to include deleted users

        :type f_email: str
        :param f_email: filter for user emails. The filter will be active for
            non-admin users only if the Galaxy instance has the
            ``expose_user_email`` option set to ``true`` in the
            ``config/galaxy.yml`` configuration file.

        :type f_name: str
        :param f_name: filter for user names. The filter will be active for
            non-admin users only if the Galaxy instance has the
            ``expose_user_name`` option set to ``true`` in the
            ``config/galaxy.yml`` configuration file.

        :type f_any: str
        :param f_any: filter for user email or name. Each filter will be active
            for non-admin users only if the Galaxy instance has the
            corresponding ``expose_user_*`` option set to ``true`` in the
            ``config/galaxy.yml`` configuration file.

        :rtype: list
        :return: a list of dicts with user details.
                 For example::

                   [{'email': 'a_user@example.org',
                     'id': 'dda47097d9189f15',
                     'url': '/api/users/dda47097d9189f15'}]

        """
        params: dict[str, Any] = {}
        if f_email:
            params["f_email"] = f_email
        if f_name:
            params["f_name"] = f_name
        if f_any:
            params["f_any"] = f_any
        return self._get(deleted=deleted, params=params)

    def show_user(self, user_id: str, deleted: bool = False) -> dict[str, Any]:
        """
        Display information about a user.

        :type user_id: str
        :param user_id: encoded user ID

        :type deleted: bool
        :param deleted: whether to return results for a deleted user

        :rtype: dict
        :return: a dictionary containing information about the user
        """
        return self._get(id=user_id, deleted=deleted)

    def create_remote_user(self, user_email: str) -> dict[str, Any]:
        """
        Create a new Galaxy remote user.

        .. note::
          This method works only if the Galaxy instance has the
          ``allow_user_creation`` and ``use_remote_user`` options set to
          ``true`` in the ``config/galaxy.yml`` configuration file. Also
          note that setting ``use_remote_user`` will require an upstream
          authentication proxy server; however, if you do not have one, access
          to Galaxy via a browser will not be possible.

        :type user_email: str
        :param user_email: email of the user to be created

        :rtype: dict
        :return: a dictionary containing information about the created user
        """
        payload = {
            "remote_user_email": user_email,
        }
        return self._post(payload)

    def create_local_user(self, username: str, user_email: str, password: str) -> dict[str, Any]:
        """
        Create a new Galaxy local user.

        .. note::
          This method works only if the Galaxy instance has the
          ``allow_user_creation`` option set to ``true`` and
          ``use_remote_user`` option set to ``false`` in the
          ``config/galaxy.yml`` configuration file.

        :type username: str
        :param username: username of the user to be created

        :type user_email: str
        :param user_email: email of the user to be created

        :type password: str
        :param password: password of the user to be created

        :rtype: dict
        :return: a dictionary containing information about the created user
        """
        payload = {
            "username": username,
            "email": user_email,
            "password": password,
        }
        return self._post(payload)

    def get_current_user(self) -> dict[str, Any]:
        """
        Display information about the user associated with this Galaxy
        connection.

        :rtype: dict
        :return: a dictionary containing information about the current user
        """
        url = self._make_url() + "/current"
        return self._get(url=url)

    def create_user_apikey(self, user_id: str) -> str:
        """
        Create a new API key for a given user.

        :type user_id: str
        :param user_id: encoded user ID

        :rtype: str
        :return: the API key for the user
        """
        url = self._make_url(user_id) + "/api_key"
        payload = {
            "user_id": user_id,
        }
        return self._post(payload, url=url)

    def delete_user(self, user_id: str, purge: bool = False) -> dict[str, Any]:
        """
        Delete a user.

        .. note::
          This method works only if the Galaxy instance has the
          ``allow_user_deletion`` option set to ``true`` in the
          ``config/galaxy.yml`` configuration file.

        :type user_id: str
        :param user_id: encoded user ID

        :type purge: bool
        :param purge: if ``True``, also purge (permanently delete) the history

        :rtype: dict
        :return: a dictionary containing information about the deleted user
        """
        params = {}
        if purge is True:
            params["purge"] = purge
        return self._delete(id=user_id, params=params)

    def get_user_apikey(self, user_id: str) -> str:
        """
        Get the current API key for a given user.

        :type user_id: str
        :param user_id: encoded user ID

        :rtype: str
        :return: the API key for the user, or 'Not available.' if it doesn't
          exist yet.
        """
        try:
            url = self._make_url(user_id) + "/api_key/detailed"
            return self._get(url=url)["key"]
        except ConnectionError as e:
            if e.status_code == 204:
                return "Not available."
            elif e.status_code != 404:
                raise
            # Galaxy 22.05 or earlier
            url = self._make_url(user_id) + "/api_key/inputs"
            return self._get(url=url)["inputs"][0]["value"]

    def get_or_create_user_apikey(self, user_id: str) -> str:
        """
        Get the current API key for a given user, creating one if it doesn't
        exist yet.

        :type user_id: str
        :param user_id: encoded user ID

        :rtype: str
        :return: the API key for the user

        .. note::
          This method works only on Galaxy 21.01 or later.
        """
        url = self._make_url(user_id) + "/api_key"
        return self._get(url=url)

    def update_user(self, user_id: str, user_data: Optional[dict] = None, **kwargs: Any) -> dict[str, Any]:
        """
        Update user information. You can either pass the attributes you want to
        change in the user_data dictionary, or provide them separately as
        keyword arguments.
        For attributes that cannot be expressed as keywords (e.g.
        extra_user_preferences use a `|` sign), pass them in user_data.

        :type user_id: str
        :param user_id: encoded user ID

        :type user_data: dict
        :param user_data: a dict containing the values to be updated, eg. { "username" : "newUsername", "email": "new@email" }

        :type username: str
        :param username: Replace user name with the given string

        :type email: str
        :param email: Replace user email with the given string

        :rtype: dict
        :return: details of the updated user
        """
        if user_data is None:
            user_data = {}
        user_data.update(kwargs)
        url = self._make_url(user_id) + "/information/inputs"
        return self._put(url=url, payload=user_data, id=user_id)
