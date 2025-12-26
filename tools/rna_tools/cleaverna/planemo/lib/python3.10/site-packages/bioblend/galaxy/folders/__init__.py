"""
Contains possible interactions with the Galaxy library folders
"""

import sys
from collections.abc import Iterator
from typing import (
    Any,
    Literal,
    Optional,
    overload,
    TYPE_CHECKING,
    Union,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class FoldersClient(Client):
    module = "folders"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def create_folder(self, parent_folder_id: str, name: str, description: Optional[str] = None) -> dict[str, Any]:
        """
        Create a folder.

        :type parent_folder_id: str
        :param parent_folder_id: Folder's description

        :type name: str
        :param name: name of the new folder

        :type description: str
        :param description: folder's description

        :rtype: dict
        :return: details of the updated folder
        """
        payload: dict[str, str] = {"name": name}
        if description:
            payload["description"] = description
        return self._post(payload=payload, id=parent_folder_id)

    @overload
    def show_folder(
        self,
        folder_id: str,
        contents: Literal[False] = False,
    ) -> dict[str, Any]: ...

    @overload
    def show_folder(
        self,
        folder_id: str,
        contents: Literal[True],
        limit: int = 10,
        offset: int = 0,
        include_deleted: bool = False,
    ) -> dict[str, Any]: ...

    def show_folder(
        self,
        folder_id: str,
        contents: bool = False,
        limit: int = 10,
        offset: int = 0,
        include_deleted: bool = False,
    ) -> dict[str, Any]:
        """
        Display information about a folder.

        :type folder_id: str
        :param folder_id: the folder's encoded id, prefixed by 'F'

        :type contents: bool
        :param contents: True to get the contents of the folder, rather
          than just the folder details.

        :type limit: int
        :param limit: When ``contents=True``, maximum number of items to return.

        :type offset: int
        :param contents: When ``contents=True``, number of items to skip. Return
          contents starting from item offset+1.

        :type include_deleted: bool
        :param include_deleted: When ``contents=True``, whether to include
          deleted items.

        :rtype: dict
        :return: dictionary including details of the folder.
          For contents=False the dict contains infos on the folder.
          For contents=True the dict contains the keys "metadata" (a dict with
          infos on the folder) and "folder_contents" (a list of dicts with info
          on the childs).

        Notes: For iterating over folder contents there is also contents_iter.
        """
        params = {
            "limit": limit,
            "offset": offset,
            "include_deleted": include_deleted,
        }
        return self._get(id=folder_id, contents=contents, params=params)

    def contents_iter(
        self,
        folder_id: str,
        batch_size: int = 10,
        include_deleted: bool = False,
    ) -> Iterator[dict[str, Any]]:
        """
        Iterate over folder contents.

        :type folder_id: str
        :param folder_id: the folder's encoded id, prefixed by 'F'

        :type batch_size: int
        :param batch_size: Batch size to be used internally.

        :type include_deleted: bool
        :param include_deleted: Whether to include deleted items.
        """
        total_rows = sys.maxsize
        params = {
            "limit": batch_size,
            "offset": 0,
            "include_deleted": include_deleted,
        }

        while params["offset"] <= total_rows:
            chunk = self._get(id=folder_id, contents=True, params=params)
            total_rows = chunk["metadata"]["total_rows"]
            yield from chunk["folder_contents"]
            params["offset"] += batch_size

    def delete_folder(self, folder_id: str, undelete: bool = False) -> dict[str, Any]:
        """
        Marks the folder with the given ``id`` as `deleted` (or removes the
        `deleted` mark if the `undelete` param is True).

        :type folder_id: str
        :param folder_id: the folder's encoded id, prefixed by 'F'

        :type undelete: bool
        :param undelete: If set to True, the folder will be undeleted
                         (i.e. the `deleted` mark will be removed)

        :return: detailed folder information
        :rtype: dict
        """
        payload = {"undelete": undelete}
        return self._delete(payload=payload, id=folder_id)

    def update_folder(self, folder_id: str, name: str, description: Optional[str] = None) -> dict[str, Any]:
        """
        Update folder information.

        :type folder_id: str
        :param folder_id: the folder's encoded id, prefixed by 'F'

        :type name: str
        :param name: name of the new folder

        :type description: str
        :param description: folder's description

        :rtype: dict
        :return: details of the updated folder
        """
        payload = {"name": name}
        if description:
            payload["description"] = description
        return self._put(payload=payload, id=folder_id)

    def get_permissions(self, folder_id: str, scope: Literal["current", "available"] = "current") -> dict[str, Any]:
        """
        Get the permissions of a folder.

        :type folder_id: str
        :param folder_id: the folder's encoded id, prefixed by 'F'

        :type scope: str
        :param scope: scope of permissions, either 'current' or 'available'

        :rtype: dict
        :return: dictionary including details of the folder permissions
        """
        url = self._make_url(folder_id) + "/permissions"
        return self._get(url=url)

    def set_permissions(
        self,
        folder_id: str,
        action: Literal["set_permissions"] = "set_permissions",
        add_ids: Optional[list[str]] = None,
        manage_ids: Optional[list[str]] = None,
        modify_ids: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Set the permissions of a folder.

        :type folder_id: str
        :param folder_id: the folder's encoded id, prefixed by 'F'

        :type action: str
        :param action: action to execute, only "set_permissions" is supported.

        :type add_ids: list of str
        :param add_ids: list of role IDs which can add datasets to the folder

        :type manage_ids: list of str
        :param manage_ids: list of role IDs which can manage datasets in the folder

        :type modify_ids: list of str
        :param modify_ids: list of role IDs which can modify datasets in the folder

        :rtype: dict
        :return: dictionary including details of the folder
        """
        url = self._make_url(folder_id) + "/permissions"
        payload: dict[str, Union[str, list[str]]] = {"action": action}
        if add_ids:
            payload["add_ids[]"] = add_ids
        if manage_ids:
            payload["manage_ids[]"] = manage_ids
        if modify_ids:
            payload["modify_ids[]"] = modify_ids
        return self._post(url=url, payload=payload)
