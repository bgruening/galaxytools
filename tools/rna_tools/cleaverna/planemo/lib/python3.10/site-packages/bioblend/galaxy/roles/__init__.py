"""
Contains possible interactions with the Galaxy Roles
"""

from typing import (
    Any,
    Optional,
    TYPE_CHECKING,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class RolesClient(Client):
    module = "roles"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_roles(self) -> list[dict[str, Any]]:
        """
        Displays a collection (list) of roles.

        :rtype: list
        :return: A list of dicts with details on individual roles.
          For example::

            [{"id": "f2db41e1fa331b3e",
              "model_class": "Role",
              "name": "Foo",
              "url": "/api/roles/f2db41e1fa331b3e"},
             {"id": "f597429621d6eb2b",
              "model_class": "Role",
              "name": "Bar",
              "url": "/api/roles/f597429621d6eb2b"}]
        """
        return self._get()

    def show_role(self, role_id: str) -> dict[str, Any]:
        """
        Display information on a single role

        :type role_id: str
        :param role_id: Encoded role ID

        :rtype: dict
        :return: Details of the given role.
          For example::

            {"description": "Private Role for Foo",
             "id": "f2db41e1fa331b3e",
             "model_class": "Role",
             "name": "Foo",
             "type": "private",
             "url": "/api/roles/f2db41e1fa331b3e"}
        """
        return self._get(id=role_id)

    def create_role(
        self,
        role_name: str,
        description: str,
        user_ids: Optional[list[str]] = None,
        group_ids: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Create a new role.

        :type role_name: str
        :param role_name: A name for the new role

        :type description: str
        :param description: Description for the new role

        :type user_ids: list
        :param user_ids: A list of encoded user IDs to add to the new role

        :type group_ids: list
        :param group_ids: A list of encoded group IDs to add to the new role

        :rtype: dict
        :return: Details of the newly created role.
          For example::

            {'description': 'desc',
             'url': '/api/roles/ebfb8f50c6abde6d',
             'model_class': 'Role',
             'type': 'admin',
             'id': 'ebfb8f50c6abde6d',
             'name': 'Foo'}

        .. versionchanged:: 0.15.0
            Changed the return value from a 1-element list to a dict.
        """
        if user_ids is None:
            user_ids = []
        if group_ids is None:
            group_ids = []
        payload = {"name": role_name, "description": description, "user_ids": user_ids, "group_ids": group_ids}
        ret = self._post(payload)
        if isinstance(ret, list):
            # Galaxy release_20.09 and earlier returned a 1-element list
            ret = ret[0]
        return ret
