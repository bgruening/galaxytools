"""
Interaction with a Tool Shed instance categories
"""

from typing import (
    Any,
    Literal,
    TYPE_CHECKING,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.toolshed import ToolShedInstance


class ToolShedCategoryClient(Client):
    module = "categories"

    def __init__(self, toolshed_instance: "ToolShedInstance") -> None:
        super().__init__(toolshed_instance)

    def get_categories(self, deleted: bool = False) -> list[dict[str, Any]]:
        """
        Returns a list of dictionaries that contain descriptions of the
        repository categories found on the given Tool Shed instance.

        :type deleted: bool
        :param deleted: whether to show deleted categories. Requires
          administrator access to the Tool Shed instance.

        :rtype: list
        :return: A list of dictionaries containing information about
          repository categories present in the Tool Shed.
          For example::

            [{'deleted': False,
              'description': 'Tools for manipulating data',
              'id': '175812cd7caaf439',
              'model_class': 'Category',
              'name': 'Text Manipulation',
              'url': '/api/categories/175812cd7caaf439'}]

        .. versionadded:: 0.5.2
        """
        return self._get(deleted=deleted)

    def show_category(self, category_id: str) -> dict[str, Any]:
        """
        Get details of a given category.

        :type category_id: str
        :param category_id: Encoded category ID

        :rtype: dict
        :return: details of the given category
        """
        return self._get(id=category_id)

    def get_repositories(
        self, category_id: str, sort_key: Literal["name", "owner"] = "name", sort_order: Literal["asc", "desc"] = "asc"
    ) -> dict[str, Any]:
        """
        Returns a dictionary of information for a repository category including
        a list of repositories belonging to the category.

        :type category_id: str
        :param category_id: Encoded category ID

        :type  sort_key: str
        :param sort_key: key for sorting. Options are 'name' or 'owner' (default 'name').

        :type  sort_order: str
        :param sort_order: ordering of sorted output. Options are 'asc' or 'desc' (default 'asc').

        :rtype: dict
        :return: A dict containing information about the category
          including a list of repository dicts.
          For example::

            {'deleted': False,
             'description': 'Tools for constructing and analyzing 3-dimensional shapes and '
                            'their properties',
             'id': '589548af7e391bcf',
             'model_class': 'Category',
             'name': 'Constructive Solid Geometry',
             'repositories': [{'create_time': '2016-08-23T18:53:23.845013',
                               'deleted': False,
                               'deprecated': False,
                               'description': 'Adds a surface field to a selected shape '
                                              'based on a given mathematical expression',
                               'homepage_url': 'https://github.com/gregvonkuster/galaxy-csg',
                               'id': 'af2ccc53697b064c',
                               'metadata': {'0:e12b55e960de': {'changeset_revision': 'e12b55e960de',
                                                               'downloadable': True,
                                                               'has_repository_dependencies': False,
                                                               'id': 'dfe022067783215f',
                                                               'includes_datatypes': False,
                                                               'includes_tool_dependencies': False,
                                                               'includes_tools': True,
                                                               'includes_tools_for_display_in_tool_panel': True,
                                                               'includes_workflows': False,
                                                               'malicious': False,
                                                               'missing_test_components': False,
                                                               'model_class': 'RepositoryMetadata',
                                                               'numeric_revision': 0,
                                                               'repository_id': 'af2ccc53697b064c'}},
                               'model_class': 'Repository',
                               'name': 'icqsol_add_surface_field_from_expression',
                               'owner': 'iuc',
                               'private': False,
                               'remote_repository_url': 'https://github.com/gregvonkuster/galaxy-csg',
                               'times_downloaded': 152,
                               'type': 'unrestricted',
                               'user_id': 'b563abc230aa8fd0'},
                              # ...
                              ],
             'repository_count': 11,
             'url': '/api/categories/589548af7e391bcf'}
        """

        params: dict[str, Any] = {}
        if sort_key:
            params.update({"sort_key": sort_key})
        if sort_order:
            params.update({"sort_order": sort_order})

        url = self._make_url(category_id) + "/repositories"
        return self._get(url=url, params=params)
