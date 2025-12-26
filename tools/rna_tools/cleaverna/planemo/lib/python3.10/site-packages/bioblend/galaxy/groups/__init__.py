"""
Contains possible interactions with the Galaxy Groups
"""

from typing import (
    Any,
    Optional,
    TYPE_CHECKING,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class GroupsClient(Client):
    module = "groups"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_groups(self) -> list[dict[str, Any]]:
        """
        Get all (not deleted) groups.

        :rtype: list
        :return: A list of dicts with details on individual groups.
          For example::

            [{'id': '33abac023ff186c2',
              'model_class': 'Group',
              'name': 'Listeria',
              'url': '/api/groups/33abac023ff186c2'},
             {'id': '73187219cd372cf8',
              'model_class': 'Group',
              'name': 'LPN',
              'url': '/api/groups/73187219cd372cf8'}]
        """
        return self._get()

    def show_group(self, group_id: str) -> dict[str, Any]:
        """
        Get details of a given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :rtype: dict
        :return: A description of group
          For example::

            {'id': '33abac023ff186c2',
             'model_class': 'Group',
             'name': 'Listeria',
             'roles_url': '/api/groups/33abac023ff186c2/roles',
             'url': '/api/groups/33abac023ff186c2',
             'users_url': '/api/groups/33abac023ff186c2/users'}
        """
        return self._get(id=group_id)

    def create_group(
        self, group_name: str, user_ids: Optional[list[str]] = None, role_ids: Optional[list[str]] = None
    ) -> list[dict[str, Any]]:
        """
        Create a new group.

        :type group_name: str
        :param group_name: A name for the new group

        :type user_ids: list
        :param user_ids: A list of encoded user IDs to add to the new group

        :type role_ids: list
        :param role_ids: A list of encoded role IDs to add to the new group

        :rtype: list
        :return: A (size 1) list with newly created group
          details, like::

            [{'id': '7c9636938c3e83bf',
              'model_class': 'Group',
              'name': 'My Group Name',
              'url': '/api/groups/7c9636938c3e83bf'}]
        """
        if user_ids is None:
            user_ids = []
        if role_ids is None:
            role_ids = []
        payload = {"name": group_name, "user_ids": user_ids, "role_ids": role_ids}
        return self._post(payload)

    def update_group(
        self,
        group_id: str,
        group_name: Optional[str] = None,
        user_ids: Optional[list[str]] = None,
        role_ids: Optional[list[str]] = None,
    ) -> None:
        """
        Update a group.

        :type group_id: str
        :param group_id: Encoded group ID

        :type group_name: str
        :param group_name: A new name for the group. If None, the group name is
          not changed.

        :type user_ids: list
        :param user_ids: New list of encoded user IDs for the group. It will
          substitute the previous list of users (with [] if not specified)

        :type role_ids: list
        :param role_ids: New list of encoded role IDs for the group. It will
          substitute the previous list of roles (with [] if not specified)

        :rtype: None
        :return: None
        """
        if user_ids is None:
            user_ids = []
        if role_ids is None:
            role_ids = []
        payload = {"name": group_name, "user_ids": user_ids, "role_ids": role_ids}
        return self._put(payload=payload, id=group_id)

    def get_group_users(self, group_id: str) -> list[dict[str, Any]]:
        """
        Get the list of users associated to the given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :rtype: list of dicts
        :return: List of group users' info
        """
        url = self._make_url(group_id) + "/users"
        return self._get(url=url)

    def get_group_roles(self, group_id: str) -> list[dict[str, Any]]:
        """
        Get the list of roles associated to the given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :rtype: list of dicts
        :return: List of group roles' info
        """
        url = self._make_url(group_id) + "/roles"
        return self._get(url=url)

    def add_group_user(self, group_id: str, user_id: str) -> dict[str, Any]:
        """
        Add a user to the given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :type user_id: str
        :param user_id: Encoded user ID to add to the group

        :rtype: dict
        :return: Added group user's info
        """
        url = "/".join((self._make_url(group_id), "users", user_id))
        return self._put(url=url)

    def add_group_role(self, group_id: str, role_id: str) -> dict[str, Any]:
        """
        Add a role to the given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :type role_id: str
        :param role_id: Encoded role ID to add to the group

        :rtype: dict
        :return: Added group role's info
        """
        url = "/".join((self._make_url(group_id), "roles", role_id))
        return self._put(url=url)

    def delete_group_user(self, group_id: str, user_id: str) -> dict[str, Any]:
        """
        Remove a user from the given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :type user_id: str
        :param user_id: Encoded user ID to remove from the group

        :rtype: dict
        :return: The user which was removed
        """
        url = "/".join((self._make_url(group_id), "users", user_id))
        return self._delete(url=url)

    def delete_group_role(self, group_id: str, role_id: str) -> dict[str, Any]:
        """
        Remove a role from the given group.

        :type group_id: str
        :param group_id: Encoded group ID

        :type role_id: str
        :param role_id: Encoded role ID to remove from the group

        :rtype: dict
        :return: The role which was removed
        """
        url = "/".join((self._make_url(group_id), "roles", role_id))
        return self._delete(url=url)
