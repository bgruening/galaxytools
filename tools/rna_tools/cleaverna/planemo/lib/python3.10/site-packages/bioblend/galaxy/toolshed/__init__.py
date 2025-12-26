"""
Interaction with a Galaxy Tool Shed.
"""

from typing import (
    Any,
    Optional,
    TYPE_CHECKING,
    Union,
)

from bioblend.galaxy.client import Client

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class ToolShedClient(Client):
    module = "tool_shed_repositories"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_repositories(self) -> list[dict[str, Any]]:
        """
        Get the list of all installed Tool Shed repositories on this Galaxy instance.

        :rtype: list
        :return: a list of dictionaries containing information about
          repositories present in the Tool Shed.
          For example::

            [{'changeset_revision': '4afe13ac23b6',
              'deleted': False,
              'dist_to_shed': False,
              'error_message': '',
              'name': 'velvet_toolsuite',
              'owner': 'edward-kirton',
              'status': 'Installed'}]

        .. versionchanged:: 0.4.1
            Changed method name from ``get_tools`` to ``get_repositories`` to
            better align with the Tool Shed concepts

        .. seealso:: bioblend.galaxy.tools.get_tool_panel()
        """
        return self._get()

    def show_repository(self, toolShed_id: str) -> dict[str, Any]:
        """
        Get details of a given Tool Shed repository as it is installed on this
        Galaxy instance.

        :type toolShed_id: str
        :param toolShed_id: Encoded Tool Shed ID

        :rtype: dict
        :return: Information about the tool
          For example::

            {'changeset_revision': 'b17455fb6222',
             'ctx_rev': '8',
             'owner': 'aaron',
             'status': 'Installed',
             'url': '/api/tool_shed_repositories/82de4a4c7135b20a'}

        .. versionchanged:: 0.4.1
            Changed method name from ``show_tool`` to ``show_repository`` to
            better align with the Tool Shed concepts
        """
        return self._get(id=toolShed_id)

    def install_repository_revision(
        self,
        tool_shed_url: str,
        name: str,
        owner: str,
        changeset_revision: str,
        install_tool_dependencies: bool = False,
        install_repository_dependencies: bool = False,
        install_resolver_dependencies: bool = False,
        tool_panel_section_id: Optional[str] = None,
        new_tool_panel_section_label: Optional[str] = None,
    ) -> Union[list[dict[str, Any]], dict[str, str]]:
        """
        Install a specified repository revision from a specified Tool Shed into
        this Galaxy instance. This example demonstrates installation of a repository
        that contains valid tools, loading them into a section of the Galaxy tool
        panel or creating a new tool panel section.
        You can choose if tool dependencies or repository dependencies should be
        installed through the Tool Shed,
        (use ``install_tool_dependencies`` or ``install_repository_dependencies``)
        or through a resolver that supports installing dependencies
        (use ``install_resolver_dependencies``). Note that any combination of
        the three dependency resolving variables is valid.

        Installing the repository into an existing tool panel section requires
        the tool panel config file (e.g., tool_conf.xml, shed_tool_conf.xml, etc)
        to contain the given tool panel section:

            <section id="from_test_tool_shed" name="From Test Tool Shed" version="">
            </section>

        :type tool_shed_url: str
        :param tool_shed_url: URL of the Tool Shed from which the repository should
                              be installed from (e.g., ``https://testtoolshed.g2.bx.psu.edu``)

        :type name: str
        :param name: The name of the repository that should be installed

        :type owner: str
        :param owner: The name of the repository owner

        :type changeset_revision: str
        :param changeset_revision: The revision of the repository to be installed

        :type install_tool_dependencies: bool
        :param install_tool_dependencies: Whether or not to automatically handle
                                          tool dependencies (see
                                          https://galaxyproject.org/toolshed/tool-dependency-recipes/
                                          for more details)

        :type install_repository_dependencies: bool
        :param install_repository_dependencies: Whether or not to automatically
                                                handle repository dependencies
                                                (see https://galaxyproject.org/toolshed/defining-repository-dependencies/
                                                for more details)

        :type install_resolver_dependencies: bool
        :param install_resolver_dependencies: Whether or not to automatically
                                                install resolver dependencies (e.g. conda).

        :type tool_panel_section_id: str
        :param tool_panel_section_id: The ID of the Galaxy tool panel section
                                      where the tool should be insterted under.
                                      Note that you should specify either this
                                      parameter or the ``new_tool_panel_section_label``.
                                      If both are specified, this one will take
                                      precedence.

        :type new_tool_panel_section_label: str
        :param new_tool_panel_section_label: The name of a Galaxy tool panel section
                                             that should be created and the repository
                                             installed into.
        """
        payload: dict[str, Any] = {}
        payload["tool_shed_url"] = tool_shed_url
        payload["name"] = name
        payload["owner"] = owner
        payload["changeset_revision"] = changeset_revision
        payload["install_tool_dependencies"] = install_tool_dependencies
        payload["install_repository_dependencies"] = install_repository_dependencies
        payload["install_resolver_dependencies"] = install_resolver_dependencies
        if tool_panel_section_id:
            payload["tool_panel_section_id"] = tool_panel_section_id
        elif new_tool_panel_section_label:
            payload["new_tool_panel_section_label"] = new_tool_panel_section_label

        url = self._make_url() + "/new/install_repository_revision"
        return self._post(url=url, payload=payload)

    def uninstall_repository_revision(
        self, name: str, owner: str, changeset_revision: str, tool_shed_url: str, remove_from_disk: bool = True
    ) -> dict[str, Any]:
        """
        Uninstalls a specified repository revision from this Galaxy instance.

        :type name: str
        :param name: The name of the repository

        :type owner: str
        :param owner: The owner of the repository

        :type changeset_revision: str
        :param changeset_revision: The revision of the repository to uninstall

        :type tool_shed_url: str
        :param tool_shed_url: URL of the Tool Shed from which the repository was
          installed from (e.g., ``https://testtoolshed.g2.bx.psu.edu``)

        :type remove_from_disk: bool
        :param remove_from_disk: whether to also remove the repository from disk
          (the default) or only deactivate it

        :rtype: dict
        :return: If successful, a dictionary with a message noting the removal
        """
        payload: dict[str, Any] = {
            "tool_shed_url": tool_shed_url,
            "name": name,
            "owner": owner,
            "changeset_revision": changeset_revision,
            "remove_from_disk": remove_from_disk,
        }
        return self._delete(params=payload)
