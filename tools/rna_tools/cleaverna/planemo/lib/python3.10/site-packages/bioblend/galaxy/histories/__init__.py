"""
Contains possible interactions with the Galaxy Histories
"""

import logging
import re
import sys
import time
import typing
import webbrowser
from re import Pattern
from typing import (
    Any,
    IO,
    Literal,
    Optional,
    overload,
    Union,
)

import bioblend
from bioblend import (
    ConnectionError,
    NotReady,
    wait_on,
)
from bioblend.galaxy.client import Client
from bioblend.galaxy.dataset_collections import CollectionDescription
from bioblend.util import attach_file

if typing.TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

log = logging.getLogger(__name__)


class HistoryClient(Client):
    gi: "GalaxyInstance"
    module = "histories"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def create_history(self, name: Optional[str] = None) -> dict[str, Any]:
        """
        Create a new history, optionally setting the ``name``.

        :type name: str
        :param name: Optional name for new history

        :rtype: dict
        :return: Dictionary containing information about newly created history
        """
        payload = {}
        if name is not None:
            payload["name"] = name
        return self._post(payload)

    def import_history(self, file_path: Optional[str] = None, url: Optional[str] = None) -> dict[str, Any]:
        """
        Import a history from an archive on disk or a URL.

        :type file_path: str
        :param file_path: Path to exported history archive on disk.

        :type url: str
        :param url: URL for an exported history archive

        :rtype: dict
        :return: Dictionary containing information about the imported history
        """
        if file_path:
            archive_file = attach_file(file_path)
            payload: dict[str, Any] = {
                "archive_source": "",
                "archive_file": archive_file,
                "archive_type": "file",
            }
        else:
            payload = {
                "archive_source": url,
                "archive_type": "url",
            }

        return self._post(payload=payload, files_attached=file_path is not None)

    def _get_histories(
        self,
        name: Optional[str] = None,
        deleted: bool = False,
        filter_user_published: Optional[bool] = None,
        get_all_published: bool = False,
        slug: Optional[str] = None,
        all: Optional[bool] = False,
        create_time_min: Optional[str] = None,
        create_time_max: Optional[str] = None,
        update_time_min: Optional[str] = None,
        update_time_max: Optional[str] = None,
        view: Optional[Literal["summary", "detailed"]] = None,
        keys: Optional[list[str]] = None,
        limit: Optional[int] = None,
        offset: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """
        Hidden method to be used by both get_histories() and get_published_histories()
        """
        assert not (filter_user_published is not None and get_all_published)

        params: dict[str, Any] = {}
        if deleted:
            params.setdefault("q", []).append("deleted")
            params.setdefault("qv", []).append(deleted)
        if filter_user_published is not None:
            params.setdefault("q", []).append("published")
            params.setdefault("qv", []).append(filter_user_published)
        if slug is not None:
            params.setdefault("q", []).append("slug")
            params.setdefault("qv", []).append(slug)
        if all:
            params["all"] = True
        if create_time_min:
            params.setdefault("q", []).append("create_time-ge")
            params.setdefault("qv", []).append(create_time_min)
        if create_time_max:
            params.setdefault("q", []).append("create_time-le")
            params.setdefault("qv", []).append(create_time_max)
        if update_time_min:
            params.setdefault("q", []).append("update_time-ge")
            params.setdefault("qv", []).append(update_time_min)
        if update_time_max:
            params.setdefault("q", []).append("update_time-le")
            params.setdefault("qv", []).append(update_time_max)
        if view:
            params["view"] = view
        if keys:
            params["keys"] = ",".join(keys)
        if limit:
            params["limit"] = limit
        if offset:
            params["offset"] = offset

        url = "/".join((self._make_url(), "published")) if get_all_published else None
        histories = self._get(url=url, params=params)

        if name is not None:
            histories = [_ for _ in histories if _["name"] == name]
        return histories

    def get_histories(
        self,
        history_id: Optional[str] = None,
        name: Optional[str] = None,
        deleted: bool = False,
        published: Optional[bool] = None,
        slug: Optional[str] = None,
        all: Optional[bool] = False,
        create_time_min: Optional[str] = None,
        create_time_max: Optional[str] = None,
        update_time_min: Optional[str] = None,
        update_time_max: Optional[str] = None,
        view: Optional[Literal["summary", "detailed"]] = None,
        keys: Optional[list[str]] = None,
        limit: Optional[int] = None,
        offset: Optional[int] = None,
    ) -> list[dict[str, Any]]:
        """
        Get all histories, or select a subset by specifying optional arguments
        for filtering (e.g. a history name).

        :type name: str
        :param name: History name to filter on.

        :type deleted: bool
        :param deleted: whether to filter for the deleted histories (``True``)
          or for the non-deleted ones (``False``)

        :type published: bool or None
        :param published: whether to filter for the published histories
          (``True``) or for the non-published ones (``False``). If not set, no
          filtering is applied. Note the filtering is only applied to the user's
          own histories; to access all histories published by any user, use the
          ``get_published_histories`` method.

        :type slug: str
        :param slug: History slug to filter on

        :type all: bool
        :param all: Whether to include histories from other users. This
          parameter works only on Galaxy 20.01 or later and can be specified
          only if the user is a Galaxy admin.

        :type create_time_min: str
        :param create_time_min: Return histories created after the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type create_time_max: str
        :param create_time_max: Return histories created before the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type update_time_min: str
        :param update_time_min: Return histories last updated after the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type update_time_max: str
        :param update_time_max: Return histories last updated before the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type view: str
        :param view: Options are 'summary' or 'detailed'. This defaults to 'summary'.
          Setting view to 'detailed' results in a larger number of fields returned.

        :type keys: List[str]
        :param keys: List of fields to return

        :type limit: int
        :param limit: Maximum number of histories to return.

        :type offset: int
        :param offset: Number of histories to skip. Return histories starting
          from item offset+1.

        :rtype: list
        :return: List of history dicts.

        .. versionchanged:: 0.17.0
            Using the deprecated ``history_id`` parameter now raises a
            ``ValueError`` exception.
        """
        if history_id is not None:
            raise ValueError(
                "The history_id parameter has been removed, use the show_history() method to view details of a history for which you know the ID.",
            )
        return self._get_histories(
            name=name,
            deleted=deleted,
            filter_user_published=published,
            get_all_published=False,
            slug=slug,
            all=all,
            create_time_min=create_time_min,
            create_time_max=create_time_max,
            update_time_min=update_time_min,
            update_time_max=update_time_max,
            view=view,
            keys=keys,
            limit=limit,
            offset=offset,
        )

    def get_published_histories(
        self,
        name: Optional[str] = None,
        deleted: bool = False,
        slug: Optional[str] = None,
        create_time_min: Optional[str] = None,
        create_time_max: Optional[str] = None,
        update_time_min: Optional[str] = None,
        update_time_max: Optional[str] = None,
    ) -> list[dict[str, Any]]:
        """
        Get all published histories (by any user), or select a subset by
        specifying optional arguments for filtering (e.g. a history name).

        :type name: str
        :param name: History name to filter on.

        :type deleted: bool
        :param deleted: whether to filter for the deleted histories (``True``)
          or for the non-deleted ones (``False``)

        :type slug: str
        :param slug: History slug to filter on

        :type create_time_min: str
        :param create_time_min: Return histories created after the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type create_time_max: str
        :param create_time_max: Return histories created before the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type update_time_min: str
        :param update_time_min: Return histories last updated after the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :type update_time_max: str
        :param update_time_max: Return histories last updated before the provided
          time and date, which should be formatted as ``YYYY-MM-DDTHH-MM-SS``.

        :rtype: list
        :return: List of history dicts.
        """
        return self._get_histories(
            name=name,
            deleted=deleted,
            filter_user_published=None,
            get_all_published=True,
            slug=slug,
            create_time_min=create_time_min,
            create_time_max=create_time_max,
            update_time_min=update_time_min,
            update_time_max=update_time_max,
        )

    @overload
    def show_history(
        self,
        history_id: str,
        contents: Literal[False] = False,
    ) -> dict[str, Any]: ...

    @overload
    def show_history(
        self,
        history_id: str,
        contents: Literal[True],
        deleted: Optional[bool] = None,
        visible: Optional[bool] = None,
        details: Optional[str] = None,
        types: Optional[list[str]] = None,
        keys: Optional[list[str]] = None,
    ) -> list[dict[str, Any]]: ...

    def show_history(
        self,
        history_id: str,
        contents: bool = False,
        deleted: Optional[bool] = None,
        visible: Optional[bool] = None,
        details: Optional[str] = None,
        types: Optional[list[str]] = None,
        keys: Optional[list[str]] = None,
    ) -> Union[dict[str, Any], list[dict[str, Any]]]:
        """
        Get details of a given history. By default, just get the history meta
        information.

        :type history_id: str
        :param history_id: Encoded history ID to filter on

        :type contents: bool
        :param contents: When ``True``, instead of the history details, return
          a list with info for all datasets in the given history.
          Note that inside each dataset info dict, the id which should be used
          for further requests about this history dataset is given by the value
          of the `id` (not `dataset_id`) key.

        :type deleted: bool or None
        :param deleted: When ``contents=True``, whether to filter for the
          deleted datasets (``True``) or for the non-deleted ones (``False``).
          If not set, no filtering is applied.

        :type visible: bool or None
        :param visible: When ``contents=True``, whether to filter for the
          visible datasets (``True``) or for the hidden ones (``False``). If not
          set, no filtering is applied.

        :type details: str
        :param details: When ``contents=True``, include dataset details. Set to
          'all' for the most information.

        :type types: list
        :param types: When ``contents=True``, filter for history content types.
          If set to ``['dataset']``, return only datasets. If set to
          ``['dataset_collection']``,  return only dataset collections. If not
          set, no filtering is applied.

        :type keys: List[str]
        :param keys: List of fields to return

        :rtype: dict or list of dicts
        :return: details of the given history or list of dataset info

        .. note::
            As an alternative to using the ``contents=True`` parameter, consider
            using ``gi.datasets.get_datasets(history_id=history_id)`` which offers
            more extensive functionality for filtering and ordering the results.

        """
        params: dict[str, Union[bool, list, str]] = {}
        if contents:
            if details:
                params["details"] = details
            if deleted is not None:
                params["deleted"] = deleted
            if visible is not None:
                params["visible"] = visible
            if types is not None:
                params["types"] = types
        if keys:
            params["keys"] = ",".join(keys)
        return self._get(id=history_id, contents=contents, params=params)

    def delete_dataset(self, history_id: str, dataset_id: str, purge: bool = False, wait: bool = False) -> None:
        """
        Mark corresponding dataset as deleted.

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_id: str
        :param dataset_id: Encoded dataset ID

        :type purge: bool
        :param purge: if ``True``, also purge (permanently delete) the dataset

        :param wait: Whether to wait for the dataset to be purged.

        :rtype: None
        :return: None

        .. note::
          The ``purge`` option works only if the Galaxy instance has the
          ``allow_user_dataset_purge`` option set to ``true`` in the
          ``config/galaxy.yml`` configuration file.
        """
        url = "/".join((self._make_url(history_id, contents=True), dataset_id))
        payload = {}
        if purge:
            payload["purge"] = True
        self._delete(payload=payload, url=url)
        if purge and wait:

            def check_dataset_purged() -> None:
                dataset = self.show_dataset(history_id, dataset_id)
                if not dataset["purged"]:
                    raise NotReady(f"Dataset {dataset_id} in library {history_id} is not purged")

            wait_on(check_dataset_purged)

    def delete_dataset_collection(self, history_id: str, dataset_collection_id: str) -> None:
        """
        Mark corresponding dataset collection as deleted.

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_collection_id: str
        :param dataset_collection_id: Encoded dataset collection ID

        :rtype: None
        :return: None
        """
        url = "/".join((self._make_url(history_id, contents=True), "dataset_collections", dataset_collection_id))
        self._delete(url=url)

    def show_dataset(self, history_id: str, dataset_id: str) -> dict[str, Any]:
        """
        Get details about a given history dataset.

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_id: str
        :param dataset_id: Encoded dataset ID

        :rtype: dict
        :return: Information about the dataset
        """
        url = "/".join((self._make_url(history_id, contents=True), dataset_id))
        return self._get(url=url)

    def show_dataset_collection(self, history_id: str, dataset_collection_id: str) -> dict[str, Any]:
        """
        Get details about a given history dataset collection.

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_collection_id: str
        :param dataset_collection_id: Encoded dataset collection ID

        :rtype: dict
        :return: Information about the dataset collection
        """
        url = "/".join((self._make_url(history_id, contents=True), "dataset_collections", dataset_collection_id))
        return self._get(url=url)

    def show_matching_datasets(
        self, history_id: str, name_filter: Optional[Union[str, Pattern[str]]] = None
    ) -> list[dict[str, Any]]:
        """
        Get dataset details for matching datasets within a history.

        :type history_id: str
        :param history_id: Encoded history ID

        :type name_filter: str
        :param name_filter: Only datasets whose name matches the
                            ``name_filter`` regular expression will be
                            returned; use plain strings for exact matches and
                            None to match all datasets in the history

        :rtype: list
        :return: List of dictionaries
        """
        if isinstance(name_filter, str):
            name_filter = re.compile(name_filter + "$")
        return [
            self.show_dataset(history_id, h["id"])
            for h in self.show_history(history_id, contents=True)
            if name_filter is None or name_filter.match(h["name"])
        ]

    def show_dataset_provenance(self, history_id: str, dataset_id: str, follow: bool = False) -> dict[str, Any]:
        """
        Get details related to how dataset was created (``id``, ``job_id``,
        ``tool_id``, ``stdout``, ``stderr``, ``parameters``, ``inputs``,
        etc...).

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_id: str
        :param dataset_id: Encoded dataset ID

        :type follow: bool
        :param follow: If ``True``, recursively fetch dataset provenance
          information for all inputs and their inputs, etc.

        :rtype: dict
        :return: Dataset provenance information
          For example::

            {'id': '6fbd9b2274c62ebe',
             'job_id': '5471ba76f274f929',
             'parameters': {'chromInfo': '"/usr/local/galaxy/galaxy-dist/tool-data/shared/ucsc/chrom/mm9.len"',
                            'dbkey': '"mm9"',
                            'experiment_name': '"H3K4me3_TAC_MACS2"',
                            'input_chipseq_file1': {'id': '6f0a311a444290f2',
                                                    'uuid': 'null'},
                            'input_control_file1': {'id': 'c21816a91f5dc24e',
                                                    'uuid': '16f8ee5e-228f-41e2-921e-a07866edce06'},
                            'major_command': '{"gsize": "2716965481.0", "bdg": "False", "__current_case__": 0, "advanced_options": {"advanced_options_selector": "off", "__current_case__": 1}, "input_chipseq_file1": 104715, "xls_to_interval": "False", "major_command_selector": "callpeak", "input_control_file1": 104721, "pq_options": {"pq_options_selector": "qvalue", "qvalue": "0.05", "__current_case__": 1}, "bw": "300", "nomodel_type": {"nomodel_type_selector": "create_model", "__current_case__": 1}}'},
             'stderr': '',
             'stdout': '',
             'tool_id': 'toolshed.g2.bx.psu.edu/repos/ziru-zhou/macs2/modencode_peakcalling_macs2/2.0.10.2',
             'uuid': '5c0c43f5-8d93-44bd-939d-305e82f213c6'}
        """
        url = "/".join((self._make_url(history_id, contents=True), dataset_id, "provenance"))
        params = {"follow": follow}
        return self._get(url=url, params=params)

    def update_history(self, history_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Update history metadata information. Some of the attributes that can be
        modified are documented below.

        :type history_id: str
        :param history_id: Encoded history ID

        :type name: str
        :param name: Replace history name with the given string

        :type annotation: str
        :param annotation: Replace history annotation with given string

        :type deleted: bool
        :param deleted: Mark or unmark history as deleted

        :type purged: bool
        :param purged: If ``True``, mark history as purged (permanently deleted).

        :type published: bool
        :param published: Mark or unmark history as published

        :type importable: bool
        :param importable: Mark or unmark history as importable

        :type tags: list
        :param tags: Replace history tags with the given list

        :rtype: dict
        :return: details of the updated history

        .. versionchanged:: 0.8.0
            Changed the return value from the status code (type int) to a dict.
        """
        return self._put(payload=kwargs, id=history_id)

    def update_dataset(self, history_id: str, dataset_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Update history dataset metadata. Some of the attributes that can be
        modified are documented below.

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_id: str
        :param dataset_id: ID of the dataset

        :type name: str
        :param name: Replace history dataset name with the given string

        :type datatype: str
        :param datatype: Replace the datatype of the history dataset with the
          given string. The string must be a valid Galaxy datatype, both the
          current and the target datatypes must allow datatype changes, and the
          dataset must not be in use as input or output of a running job
          (including uploads), otherwise an error will be raised.

        :type genome_build: str
        :param genome_build: Replace history dataset genome build (dbkey)

        :type annotation: str
        :param annotation: Replace history dataset annotation with given string

        :type deleted: bool
        :param deleted: Mark or unmark history dataset as deleted

        :type visible: bool
        :param visible: Mark or unmark history dataset as visible

        :rtype: dict
        :return: details of the updated dataset

        .. versionchanged:: 0.8.0
            Changed the return value from the status code (type int) to a dict.
        """
        url = "/".join((self._make_url(history_id, contents=True), dataset_id))
        return self._put(payload=kwargs, url=url)

    def update_dataset_collection(self, history_id: str, dataset_collection_id: str, **kwargs: Any) -> dict[str, Any]:
        """
        Update history dataset collection metadata. Some of the attributes that
        can be modified are documented below.

        :type history_id: str
        :param history_id: Encoded history ID

        :type dataset_collection_id: str
        :param dataset_collection_id: Encoded dataset_collection ID

        :type name: str
        :param name: Replace history dataset collection name with the given
          string

        :type deleted: bool
        :param deleted: Mark or unmark history dataset collection as deleted

        :type visible: bool
        :param visible: Mark or unmark history dataset collection as visible

        :rtype: dict
        :return: the updated dataset collection attributes

        .. versionchanged:: 0.8.0
            Changed the return value from the status code (type int) to a dict.
        """
        url = "/".join((self._make_url(history_id, contents=True), "dataset_collections", dataset_collection_id))
        return self._put(payload=kwargs, url=url)

    def create_history_tag(self, history_id: str, tag: str) -> dict[str, Any]:
        """
        Create history tag

        :type history_id: str
        :param history_id: Encoded history ID

        :type tag: str
        :param tag: Add tag to history

        :rtype: dict
        :return: A dictionary with information regarding the tag.
          For example::

            {'id': 'f792763bee8d277a',
             'model_class': 'HistoryTagAssociation',
             'user_tname': 'NGS_PE_RUN',
             'user_value': None}
        """
        # empty payload since we are adding the new tag using the url
        payload: dict[str, Any] = {}
        url = "/".join((self._make_url(history_id), "tags", tag))
        return self._post(payload, url=url)

    def upload_dataset_from_library(self, history_id: str, lib_dataset_id: str) -> dict[str, Any]:
        """
        Upload a dataset into the history from a library. Requires the
        library dataset ID, which can be obtained from the library
        contents.

        :type history_id: str
        :param history_id: Encoded history ID

        :type lib_dataset_id: str
        :param lib_dataset_id: Encoded library dataset ID

        :rtype: dict
        :return: Information about the newly created HDA
        """
        payload = {
            "content": lib_dataset_id,
            "source": "library",
            "from_ld_id": lib_dataset_id,  # compatibility with old API
        }
        return self._post(payload, id=history_id, contents=True)

    def create_dataset_collection(
        self,
        history_id: str,
        collection_description: Union["CollectionDescription", dict[str, Any]],
        copy_elements: bool = True,
    ) -> dict[str, Any]:
        """
        Create a new dataset collection

        :type history_id: str
        :param history_id: Encoded history ID

        :type collection_description: bioblend.galaxy.dataset_collections.CollectionDescription
        :param collection_description: a description of the dataset collection
          For example::

            {'collection_type': 'list',
             'element_identifiers': [{'id': 'f792763bee8d277a',
                                      'name': 'element 1',
                                      'src': 'hda'},
                                     {'id': 'f792763bee8d277a',
                                      'name': 'element 2',
                                      'src': 'hda'}],
             'name': 'My collection list'}

        :type copy_elements: bool
        :param copy_elements: Whether to make a copy of the elements of the
          collection being created

        :rtype: dict
        :return: Information about the new HDCA
        """
        if isinstance(collection_description, CollectionDescription):
            collection_description_dict = collection_description.to_dict()
        else:
            collection_description_dict = collection_description

        payload = {
            "name": collection_description_dict["name"],
            "type": "dataset_collection",
            "collection_type": collection_description_dict["collection_type"],
            "element_identifiers": collection_description_dict["element_identifiers"],
            "copy_elements": copy_elements,
        }
        return self._post(payload, id=history_id, contents=True)

    def delete_history(self, history_id: str, purge: bool = False) -> dict[str, Any]:
        """
        Delete a history.

        :type history_id: str
        :param history_id: Encoded history ID

        :type purge: bool
        :param purge: if ``True``, also purge (permanently delete) the history

        :rtype: dict
        :return: An error object if an error occurred or a dictionary
                 containing: ``id`` (the encoded id of the history), ``deleted`` (if the
                 history was marked as deleted), ``purged`` (if the history was
                 purged).

        .. note::
          The ``purge`` option works only if the Galaxy instance has the
          ``allow_user_dataset_purge`` option set to ``true`` in the
          ``config/galaxy.yml`` configuration file.
        """
        payload = {}
        if purge is True:
            payload["purge"] = purge
        return self._delete(payload=payload, id=history_id)

    def undelete_history(self, history_id: str) -> str:
        """
        Undelete a history

        :type history_id: str
        :param history_id: Encoded history ID

        :rtype: str
        :return: 'OK' if it was deleted
        """
        url = self._make_url(history_id, deleted=True) + "/undelete"
        return self._post(url=url)

    def get_status(self, history_id: str) -> dict[str, Any]:
        """
        Returns the state of this history

        :type history_id: str
        :param history_id: Encoded history ID

        :rtype: dict
        :return: A dict documenting the current state of the history. Has the following keys:
            'state' = This is the current state of the history, such as ok, error, new etc.
            'state_details' = Contains individual statistics for various dataset states.
            'percent_complete' = The overall number of datasets processed to completion.
        """
        history = self.show_history(history_id)

        state: dict[str, Any] = {
            "state": history["state"],
        }

        if history.get("state_details") is not None:
            state["state_details"] = history["state_details"]
            total_complete = sum(history["state_details"].values())
            if total_complete > 0:
                state["percent_complete"] = 100 * history["state_details"]["ok"] / total_complete
            else:
                state["percent_complete"] = 0
        return state

    def get_most_recently_used_history(self) -> dict[str, Any]:
        """
        Returns the current user's most recently used history (not deleted).

        :rtype: dict
        :return: History representation
        """
        url = self._make_url() + "/most_recently_used"
        return self._get(url=url)

    def export_history(
        self,
        history_id: str,
        gzip: bool = True,
        include_hidden: bool = False,
        include_deleted: bool = False,
        wait: bool = False,
        maxwait: Optional[float] = None,
    ) -> str:
        """
        Start a job to create an export archive for the given history.

        :type history_id: str
        :param history_id: history ID

        :type gzip: bool
        :param gzip: create .tar.gz archive if ``True``, else .tar

        :type include_hidden: bool
        :param include_hidden: whether to include hidden datasets
          in the export

        :type include_deleted: bool
        :param include_deleted: whether to include deleted datasets
          in the export

        :type wait: bool
        :param wait: if ``True``, block until the export is ready; else, return
          immediately

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the export to become
          ready. When set, implies that ``wait`` is ``True``.

        :rtype: str
        :return: ``jeha_id`` of the export, or empty if ``wait`` is ``False``
          and the export is not ready.
        """
        if maxwait is not None:
            assert maxwait >= 0
        else:
            if wait:
                maxwait = sys.maxsize
            else:
                maxwait = 0
        params = {
            "gzip": gzip,
            "include_hidden": include_hidden,
            "include_deleted": include_deleted,
        }
        url = f"{self._make_url(history_id)}/exports"
        time_left = maxwait
        while True:
            try:
                r = self._put(url=url, params=params)
            except ConnectionError as e:
                if e.status_code == 202:  # export is not ready
                    if time_left > 0:
                        log.info(
                            "Waiting for the export of history %s to complete. Will wait %i more s",
                            history_id,
                            time_left,
                        )
                        time.sleep(1)
                        time_left -= 1
                    else:
                        return ""
                else:
                    raise
            else:
                break
        jeha_id = r["download_url"].rsplit("/", 1)[-1]
        return jeha_id

    def download_history(
        self, history_id: str, jeha_id: str, outf: IO[bytes], chunk_size: int = bioblend.CHUNK_SIZE
    ) -> None:
        """
        Download a history export archive.  Use :meth:`export_history`
        to create an export.

        :type history_id: str
        :param history_id: history ID

        :type jeha_id: str
        :param jeha_id: jeha ID (this should be obtained via
          :meth:`export_history`)

        :type outf: file
        :param outf: output file object, open for writing in binary mode

        :type chunk_size: int
        :param chunk_size: how many bytes at a time should be read into memory

        :rtype: None
        :return: None
        """
        url = f"{self._make_url(module_id=history_id)}/exports/{jeha_id}"
        r = self.gi.make_get_request(url, stream=True)
        r.raise_for_status()
        for chunk in r.iter_content(chunk_size):
            outf.write(chunk)

    def copy_dataset(
        self, history_id: str, dataset_id: str, source: Literal["hda", "library", "library_folder"] = "hda"
    ) -> dict[str, Any]:
        """
        Copy a dataset to a history.

        :type history_id: str
        :param history_id: history ID to which the dataset should be copied

        :type dataset_id: str
        :param dataset_id: dataset ID

        :type source: str
        :param source: Source of the dataset to be copied: 'hda' (the default), 'library' or 'library_folder'

        :rtype: dict
        :return: Information about the copied dataset
        """
        return self.copy_content(history_id, dataset_id, source)

    def copy_content(
        self, history_id: str, content_id: str, source: Literal["hda", "hdca", "library", "library_folder"] = "hda"
    ) -> dict[str, Any]:
        """
        Copy existing content (e.g. a dataset) to a history.

        :type history_id: str
        :param history_id: ID of the history to which the content should be copied

        :type content_id: str
        :param content_id: ID of the content to copy

        :type source: str
        :param source: Source of the content to be copied: 'hda' (for a history
          dataset, the default), 'hdca' (for a dataset collection), 'library'
          (for a library dataset) or 'library_folder' (for all datasets in a
          library folder).

        :rtype: dict
        :return: Information about the copied content
        """

        payload = {
            "content": content_id,
            "source": source,
            "type": "dataset" if source != "hdca" else "dataset_collection",
        }

        url = self._make_url(history_id, contents=True)
        return self._post(payload=payload, url=url)

    def open_history(self, history_id: str) -> None:
        """
        Open Galaxy in a new tab of the default web browser and switch to the
        specified history.

        :type history_id: str
        :param history_id: ID of the history to switch to

        :rtype: NoneType
        :return: ``None``

        .. warning::
          After opening the specified history, all previously opened Galaxy tabs
          in the browser session will have the current history changed to this
          one, even if the interface still shows another history. Refreshing
          any such tab is recommended.
        """

        url = f"{self.gi.base_url}/history/switch_to_history?hist_id={history_id}"
        webbrowser.open_new_tab(url)

    def get_extra_files(self, history_id: str, dataset_id: str) -> list[str]:
        """
        Get extra files associated with a composite dataset, or an empty list if
        there are none.

        :type history_id: str
        :param history_id: history ID

        :type dataset_id: str
        :param dataset_id: dataset ID

        :rtype: list
        :return: List of extra files
        """
        url = "/".join((self._make_url(history_id, contents=True), dataset_id, "extra_files"))
        return self._get(url=url)
