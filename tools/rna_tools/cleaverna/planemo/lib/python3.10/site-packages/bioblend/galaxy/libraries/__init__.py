"""
Contains possible interactions with the Galaxy Data Libraries
"""

import logging
from typing import (
    Any,
    Literal,
    Optional,
    overload,
    TYPE_CHECKING,
    Union,
)

from bioblend import (
    NotReady,
    wait_on,
)
from bioblend.galaxy.client import Client
from bioblend.galaxy.datasets import TERMINAL_STATES
from bioblend.util import attach_file

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

LinkDataOnly = Literal["copy_files", "link_to_files"]

log = logging.getLogger(__name__)


class LibraryClient(Client):
    module = "libraries"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def create_library(
        self, name: str, description: Optional[str] = None, synopsis: Optional[str] = None
    ) -> dict[str, Any]:
        """
        Create a data library with the properties defined in the arguments.

        :type name: str
        :param name: Name of the new data library

        :type description: str
        :param description: Optional data library description

        :type synopsis: str
        :param synopsis: Optional data library synopsis

        :rtype: dict
        :return: Details of the created library.
          For example::

            {'id': 'f740ab636b360a70',
             'name': 'Library from bioblend',
             'url': '/api/libraries/f740ab636b360a70'}
        """
        payload = {"name": name}
        if description:
            payload["description"] = description
        if synopsis:
            payload["synopsis"] = synopsis
        return self._post(payload)

    def delete_library(self, library_id: str) -> dict[str, Any]:
        """
        Delete a data library.

        :type library_id: str
        :param library_id: Encoded data library ID identifying the library to be
          deleted

        :rtype: dict
        :return: Information about the deleted library

        .. warning::
          Deleting a data library is irreversible - all of the data from the
          library will be permanently deleted.
        """
        return self._delete(id=library_id)

    def _show_item(self, library_id: str, item_id: str) -> dict[str, Any]:
        """
        Get details about a given library item.
        """
        url = "/".join((self._make_url(library_id, contents=True), item_id))
        return self._get(url=url)

    def delete_library_dataset(self, library_id: str, dataset_id: str, purged: bool = False) -> dict[str, Any]:
        """
        Delete a library dataset in a data library.

        :type library_id: str
        :param library_id: library id where dataset is found in

        :type dataset_id: str
        :param dataset_id: id of the dataset to be deleted

        :type purged: bool
        :param purged: Indicate that the dataset should be purged (permanently
          deleted)

        :rtype: dict
        :return: A dictionary containing the dataset id and whether the dataset
          has been deleted.
          For example::

            {'deleted': True,
             'id': '60e680a037f41974'}
        """
        url = "/".join((self._make_url(library_id, contents=True), dataset_id))
        return self._delete(payload={"purged": purged}, url=url)

    def update_library_dataset(self, dataset_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Update library dataset metadata. Some of the attributes that can be
        modified are documented below.

        :type dataset_id: str
        :param dataset_id: id of the dataset to be updated

        :type name: str
        :param name: Replace library dataset name with the given string

        :type misc_info: str
        :param misc_info: Replace library dataset misc_info with given string

        :type file_ext: str
        :param file_ext: Replace library dataset extension (must exist in the Galaxy registry)

        :type genome_build: str
        :param genome_build: Replace library dataset genome build (dbkey)

        :type tags: list
        :param tags: Replace library dataset tags with the given list

        :rtype: dict
        :return: details of the updated dataset
        """
        url = "/".join((self._make_url(), "datasets", dataset_id))
        return self._patch(payload=kwargs, url=url)

    def show_dataset(self, library_id: str, dataset_id: str) -> dict[str, Any]:
        """
        Get details about a given library dataset. The required ``library_id``
        can be obtained from the datasets's library content details.

        :type library_id: str
        :param library_id: library id where dataset is found in

        :type dataset_id: str
        :param dataset_id: id of the dataset to be inspected

        :rtype: dict
        :return: A dictionary containing information about the dataset in the
          library
        """
        return self._show_item(library_id, dataset_id)

    def wait_for_dataset(
        self, library_id: str, dataset_id: str, maxwait: float = 12000, interval: float = 3
    ) -> dict[str, Any]:
        """
        Wait until the library dataset state is terminal ('ok', 'empty',
        'error', 'discarded' or 'failed_metadata').

        :type library_id: str
        :param library_id: library id where dataset is found in

        :type dataset_id: str
        :param dataset_id: id of the dataset to wait for

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the dataset state to
          become terminal. After this time, a ``TimeoutException`` will be
          raised.

        :type interval: float
        :param interval: Time (in seconds) to wait between 2 consecutive checks.

        :rtype: dict
        :return: A dictionary containing information about the dataset in the
          library
        """

        def check_and_get_library_dataset() -> dict[str, Any]:
            dataset = self.show_dataset(library_id, dataset_id)
            state = dataset["state"]
            if state in TERMINAL_STATES:
                return dataset
            raise NotReady(f"Dataset {dataset_id} in library {library_id} is in non-terminal state {state}")

        return wait_on(check_and_get_library_dataset, maxwait=maxwait, interval=interval)

    def show_folder(self, library_id: str, folder_id: str) -> dict[str, Any]:
        """
        Get details about a given folder. The required ``folder_id`` can be
        obtained from the folder's library content details.

        :type library_id: str
        :param library_id: library id to inspect folders in

        :type folder_id: str
        :param folder_id: id of the folder to be inspected

        :rtype: dict
        :return: Information about the folder
        """
        return self._show_item(library_id, folder_id)

    def _get_root_folder_id(self, library_id: str) -> str:
        """
        Find the root folder (i.e. '/') of a library.

        :type library_id: str
        :param library_id: library id to find root of
        """
        library_dict = self.show_library(library_id=library_id)
        return library_dict["root_folder_id"]

    def create_folder(
        self, library_id: str, folder_name: str, description: Optional[str] = None, base_folder_id: Optional[str] = None
    ) -> list[dict[str, Any]]:
        """
        Create a folder in a library.

        :type library_id: str
        :param library_id: library id to use

        :type folder_name: str
        :param folder_name: name of the new folder in the data library

        :type description: str
        :param description: description of the new folder in the data library

        :type base_folder_id: str
        :param base_folder_id: id of the folder where to create the new folder.
          If not provided, the root folder will be used

        :rtype: list
        :return: List with a single dictionary containing information about the new folder
        """
        # Get root folder ID if no ID was provided
        if base_folder_id is None:
            base_folder_id = self._get_root_folder_id(library_id)
        # Compose the payload
        payload = {
            "name": folder_name,
            "folder_id": base_folder_id,
            "create_type": "folder",
        }
        if description is not None:
            payload["description"] = description
        return self._post(payload, id=library_id, contents=True)

    def get_folders(
        self, library_id: str, folder_id: Optional[str] = None, name: Optional[str] = None
    ) -> list[dict[str, Any]]:
        """
        Get all the folders in a library, or select a subset by specifying a
        folder name for filtering.

        :type library_id: str
        :param library_id: library id to use

        :type name: str
        :param name: Folder name to filter on. For ``name`` specify the full
                     path of the folder starting from the library's root
                     folder, e.g. ``/subfolder/subsubfolder``.

        :rtype: list
        :return: list of dicts each containing basic information about a folder

        .. versionchanged:: 1.1.1
           Using the deprecated ``folder_id`` parameter now raises a
           ``ValueError`` exception.
        """
        if folder_id is not None:
            raise ValueError(
                "The folder_id parameter has been removed, use the show_folder() method to view details of a folder for which you know the ID."
            )
        library_contents = self.show_library(library_id=library_id, contents=True)
        if name is not None:
            folders = [_ for _ in library_contents if _["type"] == "folder" and _["name"] == name]
        else:
            folders = [_ for _ in library_contents if _["type"] == "folder"]
        return folders

    def get_libraries(
        self, library_id: Optional[str] = None, name: Optional[str] = None, deleted: Optional[bool] = False
    ) -> list[dict[str, Any]]:
        """
        Get all libraries, or select a subset by specifying optional arguments
        for filtering (e.g. a library name).

        :type name: str
        :param name: Library name to filter on.

        :type deleted: bool
        :param deleted: If ``False`` (the default), return only non-deleted
          libraries. If ``True``, return only deleted libraries. If ``None``,
          return both deleted and non-deleted libraries.

        :rtype: list
        :return: list of dicts each containing basic information about a library

        .. versionchanged:: 1.1.1
           Using the deprecated ``library_id`` parameter now raises a
           ``ValueError`` exception.
        """
        if library_id is not None:
            raise ValueError(
                "The library_id parameter has been removed, use the show_library() method to view details of a library for which you know the ID."
            )
        libraries = self._get(params={"deleted": deleted})
        if name is not None:
            libraries = [_ for _ in libraries if _["name"] == name]
        return libraries

    @overload
    def show_library(self, library_id: str, contents: Literal[False] = False) -> dict[str, Any]: ...

    @overload
    def show_library(self, library_id: str, contents: Literal[True]) -> list[dict[str, Any]]: ...

    def show_library(self, library_id: str, contents: bool = False) -> Union[dict[str, Any], list[dict[str, Any]]]:
        """
        Get information about a library.

        :type library_id: str
        :param library_id: filter for library by library id

        :type contents: bool
        :param contents: whether to get contents of the library (rather
          than just the library details)

        :return: details of the given library or list of library datasets and
          folders info
        """
        return self._get(id=library_id, contents=contents)

    def _do_upload(self, library_id: str, **kwargs: Any) -> list[dict[str, Any]]:
        """
        Set up the POST request and do the actual data upload to a data library.
        This method should not be called directly but instead refer to the
        methods specific for the desired type of data upload.
        """
        folder_id = kwargs.get("folder_id")
        if folder_id is None:
            folder_id = self._get_root_folder_id(library_id)
        files_attached = False
        # Compose the payload dict
        payload = {
            "folder_id": folder_id,
            "file_type": kwargs.get("file_type", "auto"),
            "dbkey": kwargs.get("dbkey", "?"),
            "create_type": "file",
            "tag_using_filenames": kwargs.get("tag_using_filenames", False),
            "preserve_dirs": kwargs.get("preserve_dirs", False),
        }
        if kwargs.get("roles"):
            payload["roles"] = kwargs["roles"]
        if kwargs.get("link_data_only") and kwargs["link_data_only"] != "copy_files":
            payload["link_data_only"] = "link_to_files"
        if kwargs.get("tags"):
            payload["tags"] = kwargs["tags"]
        # upload options
        if kwargs.get("file_url") is not None:
            payload["upload_option"] = "upload_file"
            payload["files_0|url_paste"] = kwargs["file_url"]
        elif kwargs.get("pasted_content") is not None:
            payload["upload_option"] = "upload_file"
            payload["files_0|url_paste"] = kwargs["pasted_content"]
        elif kwargs.get("server_dir") is not None:
            payload["upload_option"] = "upload_directory"
            payload["server_dir"] = kwargs["server_dir"]
        elif kwargs.get("file_local_path") is not None:
            payload["upload_option"] = "upload_file"
            payload["files_0|file_data"] = attach_file(kwargs["file_local_path"])
            files_attached = True
        elif kwargs.get("filesystem_paths") is not None:
            payload["upload_option"] = "upload_paths"
            payload["filesystem_paths"] = kwargs["filesystem_paths"]

        try:
            return self._post(payload, id=library_id, contents=True, files_attached=files_attached)
        finally:
            if payload.get("files_0|file_data") is not None:
                payload["files_0|file_data"].close()

    def upload_file_from_url(
        self,
        library_id: str,
        file_url: str,
        folder_id: Optional[str] = None,
        file_type: str = "auto",
        dbkey: str = "?",
        tags: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]:
        """
        Upload a file to a library from a URL.

        :type library_id: str
        :param library_id: id of the library where to place the uploaded file

        :type file_url: str
        :param file_url: URL of the file to upload

        :type folder_id: str
        :param folder_id: id of the folder where to place the uploaded file.
          If not provided, the root folder will be used

        :type file_type: str
        :param file_type: Galaxy file format name

        :type dbkey: str
        :param dbkey: Dbkey

        :type tags: list
        :param tags: A list of tags to add to the datasets

        :rtype: list
        :return: List with a single dictionary containing information about the LDDA
        """
        return self._do_upload(
            library_id, file_url=file_url, folder_id=folder_id, file_type=file_type, dbkey=dbkey, tags=tags
        )

    def upload_file_contents(
        self,
        library_id: str,
        pasted_content: str,
        folder_id: Optional[str] = None,
        file_type: str = "auto",
        dbkey: str = "?",
        tags: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]:
        """
        Upload pasted_content to a data library as a new file.

        :type library_id: str
        :param library_id: id of the library where to place the uploaded file

        :type pasted_content: str
        :param pasted_content: Content to upload into the library

        :type folder_id: str
        :param folder_id: id of the folder where to place the uploaded file.
          If not provided, the root folder will be used

        :type file_type: str
        :param file_type: Galaxy file format name

        :type dbkey: str
        :param dbkey: Dbkey

        :type tags: list
        :param tags: A list of tags to add to the datasets

        :rtype: list
        :return: List with a single dictionary containing information about the LDDA
        """
        return self._do_upload(
            library_id, pasted_content=pasted_content, folder_id=folder_id, file_type=file_type, dbkey=dbkey, tags=tags
        )

    def upload_file_from_local_path(
        self,
        library_id: str,
        file_local_path: str,
        folder_id: Optional[str] = None,
        file_type: str = "auto",
        dbkey: str = "?",
        tags: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]:
        """
        Read local file contents from file_local_path and upload data to a
        library.

        :type library_id: str
        :param library_id: id of the library where to place the uploaded file

        :type file_local_path: str
        :param file_local_path: path of local file to upload

        :type folder_id: str
        :param folder_id: id of the folder where to place the uploaded file.
          If not provided, the root folder will be used

        :type file_type: str
        :param file_type: Galaxy file format name

        :type dbkey: str
        :param dbkey: Dbkey

        :type tags: list
        :param tags: A list of tags to add to the datasets

        :rtype: list
        :return: List with a single dictionary containing information about the LDDA
        """
        return self._do_upload(
            library_id,
            file_local_path=file_local_path,
            folder_id=folder_id,
            file_type=file_type,
            dbkey=dbkey,
            tags=tags,
        )

    def upload_file_from_server(
        self,
        library_id: str,
        server_dir: str,
        folder_id: Optional[str] = None,
        file_type: str = "auto",
        dbkey: str = "?",
        link_data_only: Optional[LinkDataOnly] = None,
        roles: str = "",
        preserve_dirs: bool = False,
        tag_using_filenames: bool = False,
        tags: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]:
        """
        Upload all files in the specified subdirectory of the Galaxy library
        import directory to a library.

        :type library_id: str
        :param library_id: id of the library where to place the uploaded file

        :type server_dir: str
        :param server_dir: relative path of the subdirectory of
          ``library_import_dir`` to upload. All and only the files (i.e. no
          subdirectories) contained in the specified directory will be
          uploaded

        :type folder_id: str
        :param folder_id: id of the folder where to place the uploaded files.
          If not provided, the root folder will be used

        :type file_type: str
        :param file_type: Galaxy file format name

        :type dbkey: str
        :param dbkey: Dbkey

        :type link_data_only: str
        :param link_data_only: either 'copy_files' (default) or
          'link_to_files'. Setting to 'link_to_files' symlinks instead of
          copying the files

        :type roles: str
        :param roles: ???

        :type preserve_dirs: bool
        :param preserve_dirs: Indicate whether to preserve the directory structure when importing dir

        :type tag_using_filenames: bool
        :param tag_using_filenames: Indicate whether to generate dataset tags
          from filenames.

          .. versionchanged:: 0.14.0
            Changed the default from ``True`` to ``False``.

        :type tags: list
        :param tags: A list of tags to add to the datasets

        :rtype: list
        :return: List with a single dictionary containing information about the LDDA

        .. note::
          This method works only if the Galaxy instance has the
          ``library_import_dir`` option configured in the ``config/galaxy.yml``
          configuration file.
        """
        return self._do_upload(
            library_id,
            server_dir=server_dir,
            folder_id=folder_id,
            file_type=file_type,
            dbkey=dbkey,
            link_data_only=link_data_only,
            roles=roles,
            preserve_dirs=preserve_dirs,
            tag_using_filenames=tag_using_filenames,
            tags=tags,
        )

    def upload_from_galaxy_filesystem(
        self,
        library_id: str,
        filesystem_paths: str,
        folder_id: Optional[str] = None,
        file_type: str = "auto",
        dbkey: str = "?",
        link_data_only: Optional[LinkDataOnly] = None,
        roles: str = "",
        preserve_dirs: bool = False,
        tag_using_filenames: bool = False,
        tags: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]:
        """
        Upload a set of files already present on the filesystem of the Galaxy
        server to a library.

        :type library_id: str
        :param library_id: id of the library where to place the uploaded file

        :type filesystem_paths: str
        :param filesystem_paths: file paths on the Galaxy server to upload to
          the library, one file per line

        :type folder_id: str
        :param folder_id: id of the folder where to place the uploaded files.
          If not provided, the root folder will be used

        :type file_type: str
        :param file_type: Galaxy file format name

        :type dbkey: str
        :param dbkey: Dbkey

        :type link_data_only: str
        :param link_data_only: either 'copy_files' (default) or
          'link_to_files'. Setting to 'link_to_files' symlinks instead of
          copying the files

        :type roles: str
        :param roles: ???

        :type preserve_dirs: bool
        :param preserve_dirs: Indicate whether to preserve the directory structure when importing dir

        :type tag_using_filenames: bool
        :param tag_using_filenames: Indicate whether to generate dataset tags
          from filenames.

          .. versionchanged:: 0.14.0
            Changed the default from ``True`` to ``False``.

        :type tags: list
        :param tags: A list of tags to add to the datasets

        :rtype: list
        :return: List of dictionaries containing information about each uploaded
          LDDA.

        .. note::
          This method works only if the Galaxy instance has the
          ``allow_path_paste`` option set to ``true`` in the
          ``config/galaxy.yml`` configuration file.
        """
        return self._do_upload(
            library_id,
            filesystem_paths=filesystem_paths,
            folder_id=folder_id,
            file_type=file_type,
            dbkey=dbkey,
            link_data_only=link_data_only,
            roles=roles,
            preserve_dirs=preserve_dirs,
            tag_using_filenames=tag_using_filenames,
            tags=tags,
        )

    def copy_from_dataset(
        self, library_id: str, dataset_id: str, folder_id: Optional[str] = None, message: str = ""
    ) -> dict[str, Any]:
        """
        Copy a Galaxy dataset into a library.

        :type library_id: str
        :param library_id: id of the library where to place the uploaded file

        :type dataset_id: str
        :param dataset_id: id of the dataset to copy from

        :type folder_id: str
        :param folder_id: id of the folder where to place the uploaded files.
          If not provided, the root folder will be used

        :type message: str
        :param message: message for copying action

        :rtype: dict
        :return: LDDA information
        """
        if folder_id is None:
            folder_id = self._get_root_folder_id(library_id)
        payload = {
            "folder_id": folder_id,
            "create_type": "file",
            "from_hda_id": dataset_id,
            "ldda_message": message,
        }
        return self._post(payload, id=library_id, contents=True)

    def get_library_permissions(self, library_id: str) -> dict[str, Any]:
        """
        Get the permissions for a library.

        :type library_id: str
        :param library_id: id of the library

        :rtype: dict
        :return: dictionary with all applicable permissions' values
        """
        url = self._make_url(library_id) + "/permissions"
        return self._get(url=url)

    def get_dataset_permissions(self, dataset_id: str) -> dict[str, Any]:
        """
        Get the permissions for a dataset.

        :type dataset_id: str
        :param dataset_id: id of the dataset

        :rtype: dict
        :return: dictionary with all applicable permissions' values
        """
        url = "/".join((self._make_url(), "datasets", dataset_id, "permissions"))
        return self._get(url=url)

    def set_library_permissions(
        self,
        library_id: str,
        access_in: Optional[list[str]] = None,
        modify_in: Optional[list[str]] = None,
        add_in: Optional[list[str]] = None,
        manage_in: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Set the permissions for a library. Note: it will override all security
        for this library even if you leave out a permission type.

        :type library_id: str
        :param library_id: id of the library

        :type access_in: list
        :param access_in: list of role ids

        :type modify_in: list
        :param modify_in: list of role ids

        :type add_in: list
        :param add_in: list of role ids

        :type manage_in: list
        :param manage_in: list of role ids

        :rtype: dict
        :return: General information about the library
        """
        payload: dict[str, list[str]] = {}
        if access_in:
            payload["LIBRARY_ACCESS_in"] = access_in
        if modify_in:
            payload["LIBRARY_MODIFY_in"] = modify_in
        if add_in:
            payload["LIBRARY_ADD_in"] = add_in
        if manage_in:
            payload["LIBRARY_MANAGE_in"] = manage_in
        url = self._make_url(library_id) + "/permissions"
        return self._post(payload, url=url)

    def set_dataset_permissions(
        self,
        dataset_id: str,
        access_in: Optional[list[str]] = None,
        modify_in: Optional[list[str]] = None,
        manage_in: Optional[list[str]] = None,
    ) -> dict[str, Any]:
        """
        Set the permissions for a dataset. Note: it will override all security
        for this dataset even if you leave out a permission type.

        :type dataset_id: str
        :param dataset_id: id of the dataset

        :type access_in: list
        :param access_in: list of role ids

        :type modify_in: list
        :param modify_in: list of role ids

        :type manage_in: list
        :param manage_in: list of role ids

        :rtype: dict
        :return: dictionary with all applicable permissions' values
        """
        # we need here to define an action
        payload: dict[str, Any] = {
            "action": "set_permissions",
        }
        if access_in:
            payload["access_ids[]"] = access_in
        if modify_in:
            payload["modify_ids[]"] = modify_in
        if manage_in:
            payload["manage_ids[]"] = manage_in
        url = "/".join((self._make_url(), "datasets", dataset_id, "permissions"))
        return self._post(payload, url=url)
