"""
Interaction with a Tool Shed instance repositories
"""

from typing import (
    Any,
    Literal,
    Optional,
    TYPE_CHECKING,
    Union,
)

from bioblend.galaxy.client import Client
from bioblend.util import attach_file

if TYPE_CHECKING:
    from bioblend.toolshed import ToolShedInstance


class ToolShedRepositoryClient(Client):
    module = "repositories"

    def __init__(self, toolshed_instance: "ToolShedInstance") -> None:
        super().__init__(toolshed_instance)

    def get_repositories(self, name: Optional[str] = None, owner: Optional[str] = None) -> list[dict[str, Any]]:
        """
        Get all repositories in a Galaxy Tool Shed, or select a subset by
        specifying optional arguments for filtering (e.g. a repository name).

        :type name: str
        :param name: Repository name to filter on.

        :type owner: str
        :param owner: Repository owner to filter on.

        :rtype: list
        :return: Returns a list of dictionaries containing information about
          repositories present in the Tool Shed.
          For example::

            [{'category_ids': ['c1df3132f6334b0e', 'f6d7b0037d901d9b'],
              'create_time': '2020-02-09T16:24:37.098176',
              'deleted': False,
              'deprecated': False,
              'description': 'Order Contigs',
              'homepage_url': '',
              'id': '287bd69f724b99ce',
              'model_class': 'Repository',
              'name': 'best_tool_ever',
              'owner': 'billybob',
              'private': False,
              'remote_repository_url': '',
              'times_downloaded': 0,
              'type': 'unrestricted',
              'user_id': '5cefd48bc04af6d4'}]

        .. versionchanged:: 0.4.1
          Changed method name from ``get_tools`` to ``get_repositories`` to
          better align with the Tool Shed concepts.
        """
        params = {}
        if name:
            params["name"] = name
        if owner:
            params["owner"] = owner
        return self._get(params=params)

    def search_repositories(
        self,
        q: str,
        page: int = 1,
        page_size: int = 10,
    ) -> dict[str, Any]:
        """
        Search for repositories in a Galaxy Tool Shed.

        :type  q: str
        :param q: query string for searching purposes

        :type  page: int
        :param page: page requested

        :type  page_size: int
        :param page_size: page size requested

        :rtype:  dict
        :return: dictionary containing search hits as well as metadata for the
          search.
          For example::

            {'hits': [{'matched_terms': [],
                       'repository': {'approved': 'no',
                                      'categories': 'fastq manipulation',
                                      'description': 'Convert export file to fastq',
                                      'full_last_updated': '2015-01-18 09:48 AM',
                                      'homepage_url': '',
                                      'id': 'bdfa208f0cf6504e',
                                      'last_updated': 'less than a year',
                                      'long_description': 'This is a simple too to convert Solexas Export files to FASTQ files.',
                                      'name': 'export_to_fastq',
                                      'remote_repository_url': '',
                                      'repo_lineage': "['0:c9e926d9d87e', '1:38859774da87']"
                                      'repo_owner_username': 'louise',
                                      'times_downloaded': 164},
                       'score': 4.92},
                      {'matched_terms': [],
                       'repository': {'approved': 'no',
                                      'categories': 'fastq manipulation',
                                      'description': 'Convert BAM file to fastq',
                                      'full_last_updated': '2015-04-07 11:57 AM',
                                      'homepage_url': '',
                                      'id': '175812cd7caaf439',
                                      'last_updated': 'less than a month',
                                      'long_description': 'Use Picards SamToFastq to convert a BAM file to fastq. Useful for storing reads as BAM in Galaxy and converting to fastq when needed for analysis.',
                                      'name': 'bam_to_fastq',
                                      'remote_repository_url': '',
                                      'repo_lineage': "['0:a0af255e28c1', '1:2523cb0fb84c', '2:2656247b5253']"
                                      'repo_owner_username': 'brad-chapman',
                                      'times_downloaded': 138},
                       'score': 4.14}],
             'hostname': 'https://testtoolshed.g2.bx.psu.edu/',
             'page': '1',
             'page_size': '2',
             'total_results': '64'}
        """
        params = {"q": q, "page": page, "page_size": page_size}
        return self._get(params=params)

    def show_repository(
        self,
        toolShed_id: str,
    ) -> dict[str, Any]:
        """
        Display information of a repository from Tool Shed

        :type toolShed_id: str
        :param toolShed_id: Encoded Tool Shed ID

        :rtype: dict
        :return: Information about the tool.
          For example::

            {'category_ids': ['c1df3132f6334b0e', 'f6d7b0037d901d9b'],
             'create_time': '2020-02-22T20:39:15.548491',
             'deleted': False,
             'deprecated': False,
             'description': 'Order Contigs',
             'homepage_url': '',
             'id': '287bd69f724b99ce',
             'long_description': '',
             'model_class': 'Repository',
             'name': 'best_tool_ever',
             'owner': 'billybob',
             'private': False,
             'remote_repository_url': '',
             'times_downloaded': 0,
             'type': 'unrestricted',
             'user_id': '5cefd48bc04af6d4'}

        .. versionchanged:: 0.4.1
          Changed method name from ``show_tool`` to ``show_repository`` to
          better align with the Tool Shed concepts.
        """
        return self._get(id=toolShed_id)

    def get_ordered_installable_revisions(
        self,
        name: str,
        owner: str,
    ) -> list[str]:
        """
        Returns the ordered list of changeset revision hash strings that are
        associated with installable revisions. As in the changelog, the list is
        ordered oldest to newest.

        :type name: str
        :param name: the name of the repository

        :type owner: str
        :param owner: the owner of the repository

        :rtype: list
        :return: List of changeset revision hash strings from oldest to newest
        """
        url = self._make_url() + "/get_ordered_installable_revisions"
        params = {"name": name, "owner": owner}
        r = self._get(url=url, params=params)

        return r

    def get_repository_revision_install_info(
        self,
        name: str,
        owner: str,
        changeset_revision: str,
    ) -> list[dict[str, Any]]:
        """
        Return a list of dictionaries of metadata about a certain changeset
        revision for a single tool.

        :type name: str
        :param name: the name of the repository

        :type owner: str
        :param owner: the owner of the repository

        :type changeset_revision: str
        :param changeset_revision: the changeset_revision of the
          RepositoryMetadata object associated with the repository

        :rtype: List of dictionaries
        :return: Returns a list of the following dictionaries:

          #. a dictionary defining the repository
          #. a dictionary defining the repository revision (RepositoryMetadata)
          #. a dictionary including the additional information required to
             install the repository

          For example::

            [{'create_time': '2020-08-20T13:17:08.818518',
              'deleted': False,
              'deprecated': False,
              'description': 'Galaxy Freebayes Bayesian genetic variant detector tool',
              'homepage_url': '',
              'id': '491b7a3fddf9366f',
              'long_description': 'Galaxy Freebayes Bayesian genetic variant detector tool originally included in the Galaxy code distribution but migrated to the tool shed.',
              'model_class': 'Repository',
              'name': 'freebayes',
              'owner': 'devteam',
              'private': False,
              'remote_repository_url': '',
              'times_downloaded': 269,
              'type': 'unrestricted',
              'url': '/api/repositories/491b7a3fddf9366f',
              'user_id': '1de29d50c3c44272'},
             {'changeset_revision': 'd291dc763c4c',
              'do_not_test': False,
              'downloadable': True,
              'has_repository_dependencies': False,
              'id': '504be8aaa652c154',
              'includes_datatypes': False,
              'includes_tool_dependencies': True,
              'includes_tools': True,
              'includes_tools_for_display_in_tool_panel': True,
              'includes_workflows': False,
              'malicious': False,
              'missing_test_components': False,
              'model_class': 'RepositoryMetadata',
              'numeric_revision': 0,
              'repository_id': '491b7a3fddf9366f',
              'url': '/api/repository_revisions/504be8aaa652c154'},
              'valid_tools': [{'add_to_tool_panel': True,
                'description': '- Bayesian genetic variant detector',
                'guid': 'testtoolshed.g2.bx.psu.edu/repos/devteam/freebayes/freebayes/0.0.3',
                'id': 'freebayes',
                'name': 'FreeBayes',
                'requirements': [{'name': 'freebayes',
                  'type': 'package',
                  'version': '0.9.6_9608597d12e127c847ae03aa03440ab63992fedf'},
                {'name': 'samtools', 'type': 'package', 'version': '0.1.18'}],
                'tests': [{'inputs': [['reference_source|reference_source_selector',
                    'history'],
                  ['options_type|options_type_selector', 'basic'],
                  ['reference_source|ref_file', 'phiX.fasta'],
                  ['reference_source|input_bams_0|input_bam', 'fake_phiX_reads_1.bam']],
                  'name': 'Test-1',
                  'outputs': [['output_vcf', 'freebayes_out_1.vcf.contains']],
                  'required_files': ['fake_phiX_reads_1.bam',
                  'phiX.fasta',
                  'freebayes_out_1.vcf.contains']}],
                'tool_config': '/srv/toolshed/test/var/data/repos/000/repo_708/freebayes.xml',
                'tool_type': 'default',
                'version': '0.0.3',
                'version_string_cmd': None}]},
             {'freebayes': ['Galaxy Freebayes Bayesian genetic variant detector tool',
                            'http://testtoolshed.g2.bx.psu.edu/repos/devteam/freebayes',
                            'd291dc763c4c',
                            '9',
                            'devteam',
                            {},
                            {'freebayes/0.9.6_9608597d12e127c847ae03aa03440ab63992fedf': {'changeset_revision': 'd291dc763c4c',
                                                                                          'name': 'freebayes',
                                                                                          'repository_name': 'freebayes',
                                                                                          'repository_owner': 'devteam',
                                                                                          'type': 'package',
                                                                                          'version': '0.9.6_9608597d12e127c847ae03aa03440ab63992fedf'},
                             'samtools/0.1.18': {'changeset_revision': 'd291dc763c4c',
                                                 'name': 'samtools',
                                                 'repository_name': 'freebayes',
                                                 'repository_owner': 'devteam',
                                                 'type': 'package',
                                                 'version': '0.1.18'}}]}]
        """
        url = self._make_url() + "/get_repository_revision_install_info"
        params = {"name": name, "owner": owner, "changeset_revision": changeset_revision}
        return self._get(url=url, params=params)

    def repository_revisions(
        self,
        downloadable: Optional[bool] = None,
        malicious: Optional[bool] = None,
        missing_test_components: Optional[bool] = None,
        includes_tools: Optional[bool] = None,
    ) -> list[dict[str, Any]]:
        """
        Returns a (possibly filtered) list of dictionaries that include
        information about all repository revisions. The following parameters can
        be used to filter the list.

        :type downloadable: bool
        :param downloadable: Can the tool be downloaded

        :type malicious: bool
        :param malicious:

        :type missing_test_components: bool
        :param missing_test_components:

        :type includes_tools: bool
        :param includes_tools:

        :rtype: List of dictionaries
        :return: Returns a (possibly filtered) list of dictionaries that include
          information about all repository revisions.
          For example::

            [{'changeset_revision': '6e26c5a48e9a',
              'downloadable': True,
              'has_repository_dependencies': False,
              'id': '92250afff777a169',
              'includes_datatypes': False,
              'includes_tool_dependencies': False,
              'includes_tools': True,
              'includes_tools_for_display_in_tool_panel': True,
              'includes_workflows': False,
              'malicious': False,
              'missing_test_components': False,
              'model_class': 'RepositoryMetadata',
              'numeric_revision': None,
              'repository_id': '78f2604ff5e65707',
              'url': '/api/repository_revisions/92250afff777a169'},
             {'changeset_revision': '15a54fa11ad7',
              'downloadable': True,
              'has_repository_dependencies': False,
              'id': 'd3823c748ae2205d',
              'includes_datatypes': False,
              'includes_tool_dependencies': False,
              'includes_tools': True,
              'includes_tools_for_display_in_tool_panel': True,
              'includes_workflows': False,
              'malicious': False,
              'missing_test_components': False,
              'model_class': 'RepositoryMetadata',
              'numeric_revision': None,
              'repository_id': 'f9662009da7bfce0',
              'url': '/api/repository_revisions/d3823c748ae2205d'}]
        """
        # Not using '_make_url' or '_get' to create url since the module id used
        # to create url is not the same as needed for this method
        url = self.gi.url + "/repository_revisions"
        params = {}
        if downloadable is not None:
            params["downloadable"] = downloadable
        if malicious is not None:
            params["malicious"] = malicious
        if missing_test_components is not None:
            params["missing_test_components"] = missing_test_components
        if includes_tools is not None:
            params["includes_tools"] = includes_tools
        return self._get(url=url, params=params)

    def show_repository_revision(
        self,
        metadata_id: str,
    ) -> dict[str, Any]:
        """
        Returns a dictionary that includes information about a specified
        repository revision.

        :type metadata_id: str
        :param metadata_id: Encoded repository metadata ID

        :rtype: dict
        :return: Returns a dictionary that includes information about a
          specified repository revision.
          For example::

            {'changeset_revision': '7602de1e7f32',
             'downloadable': True,
             'has_repository_dependencies': False,
             'id': '504be8aaa652c154',
             'includes_datatypes': False,
             'includes_tool_dependencies': False,
             'includes_tools': True,
             'includes_tools_for_display_in_tool_panel': True,
             'includes_workflows': False,
             'malicious': False,
             'missing_test_components': True,
             'model_class': 'RepositoryMetadata',
             'numeric_revision': None,
             'repository_dependencies': [],
             'repository_id': '491b7a3fddf9366f',
             'url': '/api/repository_revisions/504be8aaa652c154'}
        """
        # Not using '_make_url' or '_get' to create url since the module id used
        # to create url is not the same as needed for this method
        # since metadata_id has to be defined, easy to create the url here
        url = "/".join((self.gi.url, "repository_revisions", metadata_id))
        return self._get(url=url)

    def update_repository(self, id: str, tar_ball_path: str, commit_message: Optional[str] = None) -> dict[str, Any]:
        """
        Update the contents of a Tool Shed repository with specified tar ball.

        :type id: str
        :param id: Encoded repository ID

        :type tar_ball_path: str
        :param tar_ball_path: Path to file containing tar ball to upload.

        :type commit_message: str
        :param commit_message: Commit message used for the underlying Mercurial
          repository backing Tool Shed repository.

        :rtype: dict
        :return: Returns a dictionary that includes repository content warnings.
          Most valid uploads will result in no such warning and an exception
          will be raised generally if there are problems.
          For example a successful upload will look like::

            {'content_alert': '',
             'message': ''}

        .. versionadded:: 0.5.2
        """
        url = self._make_url(id) + "/changeset_revision"
        payload: dict[str, Any] = {"file": attach_file(tar_ball_path)}
        if commit_message is not None:
            payload["commit_message"] = commit_message
        try:
            return self._post(payload=payload, files_attached=True, url=url)
        finally:
            payload["file"].close()

    def create_repository(
        self,
        name: str,
        synopsis: str,
        description: Optional[str] = None,
        type: Literal["unrestricted", "repository_suite_definition", "tool_dependency_definition"] = "unrestricted",
        remote_repository_url: Optional[str] = None,
        homepage_url: Optional[str] = None,
        category_ids: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Create a new repository in a Tool Shed.

        :type name: str
        :param name: Name of the repository

        :type synopsis: str
        :param synopsis: Synopsis of the repository

        :type description: str
        :param description: Optional description of the repository

        :type type: str
        :param type: type of the repository. One of "unrestricted",
          "repository_suite_definition", or "tool_dependency_definition"

        :type remote_repository_url: str
        :param remote_repository_url: Remote URL (e.g. GitHub/Bitbucket
          repository)

        :type homepage_url: str
        :param homepage_url: Upstream's homepage for the project

        :type category_ids: list
        :param category_ids: List of encoded category IDs

        :rtype: dict
        :return: a dictionary containing information about the new repository.
          For example::

            {"deleted": false,
             "deprecated": false,
             "description": "new_synopsis",
             "homepage_url": "https://github.com/galaxyproject/",
             "id": "8cf91205f2f737f4",
             "long_description": "this is some repository",
             "model_class": "Repository",
             "name": "new_repo_17",
             "owner": "qqqqqq",
             "private": false,
             "remote_repository_url": "https://github.com/galaxyproject/tools-devteam",
             "times_downloaded": 0,
             "type": "unrestricted",
             "user_id": "adb5f5c93f827949"}
        """
        payload: dict[str, Union[str, list[str]]] = {
            "name": name,
            "synopsis": synopsis,
        }
        if description is not None:
            payload["description"] = description
        if description is not None:
            payload["description"] = description
        if type is not None:
            payload["type"] = type
        if remote_repository_url is not None:
            payload["remote_repository_url"] = remote_repository_url
        if homepage_url is not None:
            payload["homepage_url"] = homepage_url
        if category_ids is not None:
            payload["category_ids[]"] = category_ids
        return self._post(payload)

    def update_repository_metadata(
        self,
        toolShed_id: str,
        name: Optional[str] = None,
        synopsis: Optional[str] = None,
        description: Optional[str] = None,
        remote_repository_url: Optional[str] = None,
        homepage_url: Optional[str] = None,
        category_ids: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Update metadata of a Tool Shed repository.

        :type toolShed_id: str
        :param name: ID of the repository to update

        :type name: str
        :param name: New name of the repository

        :type synopsis: str
        :param synopsis: New synopsis of the repository

        :type description: str
        :param description: New description of the repository

        :type remote_repository_url: str
        :param remote_repository_url: New remote URL (e.g. GitHub/Bitbucket
          repository)

        :type homepage_url: str
        :param homepage_url: New upstream homepage for the project

        :type category_ids: list
        :param category_ids: New list of encoded category IDs

        :rtype: dict
        :return: a dictionary containing information about the updated repository.
        """
        payload: dict[str, Union[str, list[str]]] = {}
        if name:
            payload["name"] = name
        if synopsis:
            payload["synopsis"] = synopsis
        if description:
            payload["description"] = description
        if remote_repository_url:
            payload["remote_repository_url"] = remote_repository_url
        if homepage_url:
            payload["homepage_url"] = homepage_url
        if category_ids:
            payload["category_ids"] = category_ids
        return self._put(id=toolShed_id, payload=payload)
