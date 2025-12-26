"""
Contains possible interaction dealing with Galaxy tools.
"""

from os.path import basename
from typing import (
    Any,
    Literal,
    Optional,
    TYPE_CHECKING,
    Union,
)

from bioblend.galaxy.client import Client
from bioblend.galaxyclient import UPLOAD_CHUNK_SIZE
from bioblend.util import attach_file
from .inputs import InputsBuilder

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance


class ToolClient(Client):
    gi: "GalaxyInstance"
    module = "tools"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def get_tools(
        self, tool_id: Optional[str] = None, name: Optional[str] = None, trackster: Optional[bool] = None
    ) -> list[dict[str, Any]]:
        """
        Get all tools, or select a subset by specifying optional arguments for
        filtering (e.g. a tool name).

        :type name: str
        :param name: Tool name to filter on.

        :type trackster: bool
        :param trackster: whether to return only tools that are compatible with
          Trackster

        :rtype: list
        :return: List of tool descriptions.

        .. seealso:: bioblend.galaxy.toolshed.get_repositories()

        .. versionchanged:: 1.1.1
           Using the deprecated ``tool_id`` parameter now raises a
           ``ValueError`` exception.
        """
        if tool_id is not None:
            raise ValueError(
                "The tool_id parameter has been removed, use the show_tool() method to view details of a tool for which you know the ID."
            )
        tools = self._raw_get_tool(in_panel=False, trackster=trackster)
        if name is not None:
            tools = [_ for _ in tools if _["name"] == name]
        return tools

    def get_tool_panel(self) -> list[dict[str, Any]]:
        """
        Get a list of available tool elements in Galaxy's configured toolbox.

        :rtype: list
        :return: List containing tools (if not in sections) or tool sections
                 with nested tool descriptions.

        .. seealso:: bioblend.galaxy.toolshed.get_repositories()
        """
        return self._raw_get_tool(in_panel=True)

    def _raw_get_tool(self, in_panel: Optional[bool] = None, trackster: Optional[bool] = None) -> list[dict[str, Any]]:
        params = {
            "in_panel": in_panel,
            "trackster": trackster,
        }
        return self._get(params=params)

    def requirements(self, tool_id: str) -> list[dict[str, Any]]:
        """
        Return the resolver status for a specific tool.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :rtype: list
        :return: List containing a resolver status dict for each tool
          requirement. For example::

            [{'cacheable': False,
              'dependency_resolver': {'auto_init': True,
                                      'auto_install': False,
                                      'can_uninstall_dependencies': True,
                                      'ensure_channels': 'iuc,conda-forge,bioconda,defaults',
                                      'model_class': 'CondaDependencyResolver',
                                      'prefix': '/mnt/galaxy/tool_dependencies/_conda',
                                      'resolver_type': 'conda',
                                      'resolves_simple_dependencies': True,
                                      'use_local': False,
                                      'versionless': False},
              'dependency_type': 'conda',
              'environment_path': '/mnt/galaxy/tool_dependencies/_conda/envs/__blast@2.10.1',
              'exact': True,
              'model_class': 'MergedCondaDependency',
              'name': 'blast',
              'version': '2.10.1'}]

        .. note::
          This method works only if the user is a Galaxy admin.
        """
        url = self._make_url(tool_id) + "/requirements"
        return self._get(url=url)

    def reload(self, tool_id: str) -> dict:
        """
        Reload the specified tool in the toolbox.

        Any changes that have been made to the wrapper since the tool was last
        reloaded will take effect.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :rtype: dict
        :return: dict containing the id, name, and version of the reloaded tool.
          For example::

            {'message': {'id': 'toolshed.g2.bx.psu.edu/repos/lparsons/cutadapt/cutadapt/3.4+galaxy1',
                         'name': 'Cutadapt',
                         'version': '3.4+galaxy1'}}

        .. note::
          This method works only if the user is a Galaxy admin.
        """
        url = self._make_url(tool_id) + "/reload"
        return self._put(url=url)

    def get_citations(self, tool_id: str) -> list[dict]:
        """
        Get BibTeX citations for a given tool ID.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :rtype: list of dicts
        :return: list containing the citations
        """
        url = self._make_url(tool_id) + "/citations"
        return self._get(url=url)

    def install_dependencies(self, tool_id: str) -> list[dict[str, Any]]:
        """
        Install dependencies for a given tool via a resolver.
        This works only for Conda currently.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :rtype: list of dicts
        :return: List of tool requirement status dictionaries

        .. note::
          This method works only if the user is a Galaxy admin.
        """
        url = self._make_url(tool_id) + "/install_dependencies"
        return self._post(url=url)

    def uninstall_dependencies(self, tool_id: str) -> dict:
        """
        Uninstall dependencies for a given tool via a resolver.
        This works only for Conda currently.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :rtype: dict
        :return: Tool requirement status

        .. note::
          This method works only if the user is a Galaxy admin.
        """
        url = self._make_url(tool_id) + "/dependencies"
        return self._delete(url=url)

    def show_tool(self, tool_id: str, io_details: bool = False, link_details: bool = False) -> dict[str, Any]:
        """
        Get details of a given tool.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :type io_details: bool
        :param io_details: whether to get also input and output details

        :type link_details: bool
        :param link_details: whether to get also link details

        :rtype: dict
        :return: Information about the tool's interface
        """
        params = {
            "io_details": io_details,
            "link_details": link_details,
        }
        return self._get(id=tool_id, params=params)

    def get_tool_tests(self, tool_id: str, tool_version: Optional[str] = None) -> list[dict[str, Any]]:
        """
        Fetch the test case definitions configured for a tool.

        :type tool_id: str
        :param tool_id: id of the requested tool

        :type tool_version: str
        :param tool_version: optional tool version selector. Use ``'*'`` to retrieve
          the tests for all available tool versions.

        :rtype: list
        :return: List of test case definitions, e.g.::

            [
                {
                    "name": "Test-1",
                    "tool_id": "random_lines1",
                    "tool_version": "2.0.2",
                    "inputs": {"num_lines": ["1"], "input": ["1.bed"]},
                    "outputs": [{"name": "out_file1", "value": "1.bed"}],
                    "required_files": [["1.bed", {"value": "1.bed"}]],
                },
                ...
            ]
        """
        params: dict[str, str] = {}
        if tool_version is not None:
            params["tool_version"] = tool_version
        url = self._make_url(tool_id) + "/test_data"
        return self._get(url=url, params=params)

    def build(
        self,
        tool_id: str,
        inputs: Optional[dict[str, Any]] = None,
        tool_version: Optional[str] = None,
        history_id: Optional[str] = None,
    ) -> dict[str, Any]:
        """
        This method returns the tool model, which includes an updated input parameter array for the given tool,
        based on user-defined "inputs".

        :type inputs: dict
        :param inputs: (optional) inputs for the payload.
          For example::

            {
                "num_lines": "1",
                "input": {
                    "values": [
                        {
                            "src": "hda",
                            "id": "4d366c1196c36d18"
                        }
                    ]
                },
                "seed_source|seed_source_selector": "no_seed",
            }

        :type tool_id: str
        :param tool_id: id of the requested tool

        :type history_id: str
        :param history_id: id of the requested history

        :type tool_version: str
        :param tool_version: version of the requested tool

        :rtype: dict
        :return: Returns a tool model including dynamic parameters and updated values, repeats block etc.
          For example::

            {
                "model_class": "Tool",
                "id": "random_lines1",
                "name": "Select random lines",
                "version": "2.0.2",
                "description": "from a file",
                "labels": [],
                "edam_operations": [],
                "edam_topics": [],
                "hidden": "",
                "is_workflow_compatible": True,
                "xrefs": [],
                "config_file": "/Users/joshij/galaxy/tools/filters/randomlines.xml",
                "panel_section_id": "textutil",
                "panel_section_name": "Text Manipulation",
                "form_style": "regular",
                "inputs": [
                    {
                        "model_class": "IntegerToolParameter",
                        "name": "num_lines",
                        "argument": None,
                        "type": "integer",
                        "label": "Randomly select",
                        "help": "lines",
                        "refresh_on_change": False,
                        "min": None,
                        "max": None,
                        "optional": False,
                        "hidden": False,
                        "is_dynamic": False,
                        "value": "1",
                        "area": False,
                        "datalist": [],
                        "default_value": "1",
                        "text_value": "1",
                    },
                ],
                "help": 'This tool selects N random lines from a file, with no repeats, and preserving ordering.',
                "citations": False,
                "sharable_url": None,
                "message": "",
                "warnings": "",
                "versions": ["2.0.2"],
                "requirements": [],
                "errors": {},
                "tool_errors": None,
                "state_inputs": {
                    "num_lines": "1",
                    "input": {"values": [{"id": "4d366c1196c36d18", "src": "hda"}]},
                    "seed_source": {"seed_source_selector": "no_seed", "__current_case__": 0},
                },
                "job_id": None,
                "job_remap": None,
                "history_id": "c9468fdb6dc5c5f1",
                "display": True,
                "action": "/tool_runner/index",
                "license": None,
                "creator": None,
                "method": "post",
                "enctype": "application/x-www-form-urlencoded",
            }

        """
        params: dict[str, Union[str, dict]] = {}

        if inputs:
            params["inputs"] = inputs

        if tool_version:
            params["tool_version"] = tool_version

        if history_id:
            params["history_id"] = history_id

        url = "/".join((self.gi.url, "tools", tool_id, "build"))

        return self._post(payload=params, url=url)

    def run_tool(
        self,
        history_id: str,
        tool_id: str,
        tool_inputs: Union[InputsBuilder, dict],
        input_format: Literal["21.01", "legacy"] = "legacy",
        data_manager_mode: Optional[Literal["populate", "dry_run", "bundle"]] = None,
    ) -> dict[str, Any]:
        """
        Runs tool specified by ``tool_id`` in history indicated
        by ``history_id`` with inputs from ``dict`` ``tool_inputs``.

        :type history_id: str
        :param history_id: encoded ID of the history in which to run the tool

        :type tool_id: str
        :param tool_id: ID of the tool to be run

        :type data_manager_mode: str
        :param data_manager_mode: Possible values are 'populate', 'dry_run' and 'bundle'.

          'populate' is the default behavior for data manager tools and results in tool data table
          files being updated after the data manager job completes.

          'dry_run' will skip any processing after the data manager job completes

          'bundle' will create a data manager bundle that can be imported on other Galaxy servers.

        :type tool_inputs: dict
        :param tool_inputs: dictionary of input datasets and parameters
          for the tool (see below)

        :type input_format:  string
        :param input_format: input format for the payload. Possible values are the
          default 'legacy' (where inputs nested inside conditionals
          or repeats are identified with e.g. '<conditional_name>|<input_name>')
          or '21.01' (where inputs inside conditionals or repeats are nested elements).

        :rtype: dict
        :return: Information about outputs and job
          For example::

            {'implicit_collections': [],
             'jobs': [{'create_time': '2019-05-08T12:26:16.067372',
                       'exit_code': None,
                       'id': '7dd125b61b35d782',
                       'model_class': 'Job',
                       'state': 'new',
                       'tool_id': 'cut1',
                       'update_time': '2019-05-08T12:26:16.067389'}],
             'output_collections': [],
             'outputs': [{'create_time': '2019-05-08T12:26:15.997739',
                          'data_type': 'galaxy.datatypes.tabular.Tabular',
                          'deleted': False,
                          'file_ext': 'tabular',
                          'file_size': 0,
                          'genome_build': '?',
                          'hda_ldda': 'hda',
                          'hid': 42,
                          'history_content_type': 'dataset',
                          'history_id': 'df8fe5ddadbf3ab1',
                          'id': 'aeb65580396167f3',
                          'metadata_column_names': None,
                          'metadata_column_types': None,
                          'metadata_columns': None,
                          'metadata_comment_lines': None,
                          'metadata_data_lines': None,
                          'metadata_dbkey': '?',
                          'metadata_delimiter': '\t',
                          'misc_blurb': 'queued',
                          'misc_info': None,
                          'model_class': 'HistoryDatasetAssociation',
                          'name': 'Cut on data 1',
                          'output_name': 'out_file1',
                          'peek': None,
                          'purged': False,
                          'state': 'new',
                          'tags': [],
                          'update_time': '2019-05-08T12:26:16.069798',
                          'uuid': 'd91d10af-7546-45be-baa9-902010661466',
                          'visible': True}]}

        The ``tool_inputs`` dict should contain input datasets and parameters
        in the (largely undocumented) format used by the Galaxy API. If you are unsure
        how to construct this dict for the tool you want to run, you can obtain a
        template by executing the ``build()`` method and taking the value of
        ``state_inputs`` from its output, then modifying it as you require.
        You can also check the examples in `Galaxy's API test suite
        <https://github.com/galaxyproject/galaxy/blob/dev/lib/galaxy_test/api/test_tools.py>`_.
        """
        payload: dict[str, Union[str, dict]] = {
            "history_id": history_id,
            "tool_id": tool_id,
            "input_format": input_format,
        }

        if isinstance(tool_inputs, InputsBuilder):
            payload["inputs"] = tool_inputs.to_dict()
        else:
            payload["inputs"] = tool_inputs

        if data_manager_mode:
            payload["data_manager_mode"] = data_manager_mode

        return self._post(payload)

    def upload_file(
        self,
        path: str,
        history_id: str,
        storage: Optional[str] = None,
        metadata: Optional[dict] = None,
        chunk_size: Optional[int] = UPLOAD_CHUNK_SIZE,
        **kwargs: Any,
    ) -> dict[str, Any]:
        """
        Upload the file specified by ``path`` to the history specified by
        ``history_id``.

        :type path: str
        :param path: path of the file to upload

        :type history_id: str
        :param history_id: id of the history where to upload the file

        :type storage: str
        :param storage: Local path to store URLs resuming uploads

        :type metadata: dict
        :param metadata: Metadata to send with upload request

        :type chunk_size: int
        :param chunk_size: Number of bytes to send in each chunk

        :type file_name: str
        :param file_name: (optional) name of the new history dataset

        :type file_type: str
        :param file_type: (optional) Galaxy datatype for the new dataset, default is auto

        :type dbkey: str
        :param dbkey: (optional) genome dbkey

        :type to_posix_lines: bool
        :param to_posix_lines: if ``True`` (the default), convert universal line
          endings to POSIX line endings. Set to ``False`` when uploading a gzip,
          bz2 or zip archive containing a binary file

        :type space_to_tab: bool
        :param space_to_tab: whether to convert spaces to tabs. Default is
          ``False``. Applicable only if to_posix_lines is ``True``

        :type auto_decompress: bool
        :param auto_decompress: Automatically decompress files if the uploaded
          file is compressed and the file type is not one that supports
          compression (e.g. ``fastqsanger.gz``). Default is ``False``.

        :rtype: dict
        :return: Information about the created upload job

        .. note::
          The following parameters work only on Galaxy 22.01 or later:
          ``storage``, ``metadata``, ``chunk_size``, ``auto_decompress``.
        """
        if self.gi.config.get_version()["version_major"] >= "22.01":
            # Use the tus protocol
            uploader = self.gi.get_tus_uploader(path, storage=storage, metadata=metadata, chunk_size=chunk_size)
            uploader.upload()
            return self.post_to_fetch(path, history_id, uploader.session_id, **kwargs)
        else:
            if "file_name" not in kwargs:
                kwargs["file_name"] = basename(path)
            payload = self._upload_payload(history_id, **kwargs)
            payload["files_0|file_data"] = attach_file(path, name=kwargs["file_name"])
            try:
                return self._post(payload, files_attached=True)
            finally:
                payload["files_0|file_data"].close()

    def post_to_fetch(self, path: str, history_id: str, session_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Make a POST request to the Fetch API after performing a tus upload.

        This is called by :meth:`upload_file` after performing an upload. This
        method is useful if you want to control the tus uploader yourself (e.g.
        to report on progress)::

            uploader = gi.get_tus_uploader(path, storage=storage)
            while uploader.offset < uploader.file_size:
                uploader.upload_chunk()
                # perform other actions...
            gi.tools.post_to_fetch(path, history_id, uploader.session_id, **upload_kwargs)

        :type session_id: str
        :param session_id: Session ID returned by the tus service

        See :meth:`upload_file` for additional parameters.

        :rtype: dict
        :return: Information about the created upload job
        """
        payload = self._fetch_payload(path, history_id, session_id, **kwargs)
        url = "/".join((self.gi.url, "tools/fetch"))
        return self._post(payload, url=url)

    def upload_from_ftp(self, path: str, history_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Upload the file specified by ``path`` from the user's FTP directory to
        the history specified by ``history_id``.

        :type path: str
        :param path: path of the file in the user's FTP directory

        :type history_id: str
        :param history_id: id of the history where to upload the file

        See :meth:`upload_file` for the optional parameters.

        :rtype: dict
        :return: Information about the created upload job
        """
        payload = self._upload_payload(history_id, **kwargs)
        payload["files_0|ftp_files"] = path
        return self._post(payload)

    def paste_content(self, content: str, history_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Upload a string to a new dataset in the history specified by
        ``history_id``.

        :type content: str
        :param content: content of the new dataset to upload or a list of URLs
          (one per line) to upload

        :type history_id: str
        :param history_id: id of the history where to upload the content

        :rtype: dict
        :return: Information about the created upload job

        See :meth:`upload_file` for the optional parameters.
        """
        payload = self._upload_payload(history_id, **kwargs)
        payload["files_0|url_paste"] = content
        return self._post(payload, files_attached=False)

    put_url = paste_content

    def _upload_payload(self, history_id: str, **kwargs: Any) -> dict[str, Any]:
        tool_input: dict[str, Any] = {
            "file_type": kwargs.get("file_type", "auto"),
            "dbkey": kwargs.get("dbkey", "?"),
            "files_0|type": "upload_dataset",
        }
        if not kwargs.get("to_posix_lines", True):
            tool_input["files_0|to_posix_lines"] = False
        elif kwargs.get("space_to_tab", False):
            tool_input["files_0|space_to_tab"] = "Yes"
        if "file_name" in kwargs:
            tool_input["files_0|NAME"] = kwargs["file_name"]
        return {
            "history_id": history_id,
            "tool_id": kwargs.get("tool_id", "upload1"),
            "inputs": tool_input,
        }

    def _fetch_payload(self, path: str, history_id: str, session_id: str, **kwargs: Any) -> dict:
        file_name = kwargs.get("file_name", basename(path))
        element = {
            "src": "files",
            "ext": kwargs.get("file_type", "auto"),
            "dbkey": kwargs.get("dbkey", "?"),
            "to_posix_lines": kwargs.get("to_posix_lines", True),
            "space_to_tab": kwargs.get("space_to_tab", False),
            "name": file_name,
        }
        payload = {
            "history_id": history_id,
            "targets": [
                {
                    "destination": {"type": "hdas"},
                    "elements": [element],
                }
            ],
            "files_0|file_data": {"session_id": session_id, "name": file_name},
            "auto_decompress": kwargs.get("auto_decompress", False),
        }
        return payload
