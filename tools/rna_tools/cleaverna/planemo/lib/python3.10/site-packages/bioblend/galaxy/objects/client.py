"""
Clients for interacting with specific Galaxy entity types.

Classes in this module should not be instantiated directly, but used
via their handles in :class:`~.galaxy_instance.GalaxyInstance`.
"""

import abc
import builtins
import json
from collections.abc import Sequence
from typing import (
    Any,
    cast,
    Generic,
    Literal,
    Optional,
    overload,
    TYPE_CHECKING,
    Union,
)

import bioblend
from bioblend.galaxy.datasets import HdaLdda
from . import wrappers

if TYPE_CHECKING:
    from .galaxy_instance import GalaxyInstance


class ObjClient(abc.ABC):
    def __init__(self, obj_gi: "GalaxyInstance") -> None:
        self.obj_gi = obj_gi
        self.gi = self.obj_gi.gi
        self.log = bioblend.log

    @abc.abstractmethod
    def get(self, id_: str) -> wrappers.Wrapper:
        """
        Retrieve the object corresponding to the given id.
        """

    @abc.abstractmethod
    def get_previews(self, **kwargs: Any) -> list:
        """
        Get a list of object previews.

        Previews entity summaries provided by REST collection URIs, e.g.
        ``http://host:port/api/libraries``.  Being the most lightweight objects
        associated to the various entities, these are the ones that should be
        used to retrieve their basic info.

        :rtype: list
        :return: a list of object previews
        """

    @abc.abstractmethod
    def list(self) -> list:
        """
        Get a list of objects.

        This method first gets the entity summaries, then gets the complete
        description for each entity with an additional GET call, so may be slow.

        :rtype: list
        :return: a list of objects
        """

    def _select_id(self, id_: Optional[str] = None, name: Optional[str] = None) -> str:
        """
        Return the id that corresponds to the given id or name info.
        """
        if id_ is None and name is None:
            raise ValueError("Neither id nor name provided")
        if id_ is not None and name is not None:
            raise ValueError("Both id and name provided")
        if id_ is None:
            id_list = [_.id for _ in self.get_previews(name=name)]
            if len(id_list) > 1:
                raise ValueError(f"Ambiguous name '{name}'")
            if not id_list:
                raise ValueError(f"Name '{name}' not found")
            return id_list[0]
        else:
            return id_

    def _get_dict(
        self, meth_name: str, reply: Optional[Union[dict[str, Any], builtins.list[dict[str, Any]]]]
    ) -> dict[str, Any]:
        if reply is None:
            raise RuntimeError(f"{meth_name}: no reply")
        elif isinstance(reply, dict):
            return reply
        try:
            return reply[0]
        except (TypeError, IndexError):
            raise RuntimeError(f"{meth_name}: unexpected reply: {reply!r}")


class ObjDatasetContainerClient(
    ObjClient, Generic[wrappers.DatasetContainerSubtype, wrappers.DatasetContainerPreviewSubtype]
):
    CONTAINER_TYPE: type[wrappers.DatasetContainer]
    CONTAINER_PREVIEW_TYPE: type[wrappers.DatasetContainerPreview]

    def __init__(self, obj_gi: "GalaxyInstance") -> None:
        super().__init__(obj_gi=obj_gi)
        gi_client = getattr(self.gi, self.CONTAINER_TYPE.API_MODULE)
        get_fname = f"get_{self.CONTAINER_TYPE.API_MODULE}"
        self._get_f = getattr(gi_client, get_fname)
        show_fname = f"show_{self.CONTAINER_TYPE.__name__.lower()}"
        self._show_f = getattr(gi_client, show_fname)

    def get_previews(
        self, name: Optional[str] = None, deleted: bool = False, **kwargs: Any
    ) -> list[wrappers.DatasetContainerPreviewSubtype]:
        dicts = self._get_f(name=name, deleted=deleted, **kwargs)
        return [
            cast(wrappers.DatasetContainerPreviewSubtype, self.CONTAINER_PREVIEW_TYPE(_, gi=self.obj_gi)) for _ in dicts
        ]

    def get(self, id_: str) -> wrappers.DatasetContainerSubtype:
        """
        Retrieve the dataset container corresponding to the given id.
        """
        cdict = self._show_f(id_)
        c_infos = self._show_f(id_, contents=True)
        if not isinstance(c_infos, Sequence):
            raise RuntimeError(f"{self._show_f.__name__}: unexpected reply: {c_infos!r}")
        c_infos = [self.CONTAINER_TYPE.CONTENT_INFO_TYPE(_) for _ in c_infos]
        return cast(wrappers.DatasetContainerSubtype, self.CONTAINER_TYPE(cdict, content_infos=c_infos, gi=self.obj_gi))


class ObjLibraryClient(ObjDatasetContainerClient[wrappers.Library, wrappers.LibraryPreview]):
    """
    Interacts with Galaxy libraries.
    """

    CONTAINER_TYPE = wrappers.Library
    CONTAINER_PREVIEW_TYPE = wrappers.LibraryPreview

    def create(self, name: str, description: Optional[str] = None, synopsis: Optional[str] = None) -> wrappers.Library:
        """
        Create a data library with the properties defined in the arguments.

        :rtype: :class:`~.wrappers.Library`
        :return: the library just created
        """
        res = self.gi.libraries.create_library(name, description, synopsis)
        lib_info = self._get_dict("create_library", res)
        return self.get(lib_info["id"])

    def list(self, name: Optional[str] = None, deleted: bool = False) -> list[wrappers.Library]:
        """
        Get libraries owned by the user of this Galaxy instance.

        :type name: str
        :param name: return only libraries with this name
        :type deleted: bool
        :param deleted: if ``True``, return libraries that have been deleted

        :rtype: list of :class:`~.wrappers.Library`
        """
        dicts = self.gi.libraries.get_libraries(name=name, deleted=deleted)
        if not deleted:
            # return Library objects only for not-deleted libraries since Galaxy
            # does not filter them out and Galaxy release_14.08 and earlier
            # crashes when trying to get a deleted library
            return [self.get(_["id"]) for _ in dicts if not _["deleted"]]
        else:
            return [self.get(_["id"]) for _ in dicts]

    def delete(self, id_: Optional[str] = None, name: Optional[str] = None) -> None:
        """
        Delete the library with the given id or name.

        Fails if multiple libraries have the specified name.

        .. warning::
          Deleting a data library is irreversible - all of the data from
          the library will be permanently deleted.
        """
        id_ = self._select_id(id_=id_, name=name)
        res = self.gi.libraries.delete_library(id_)
        if not isinstance(res, dict):
            raise RuntimeError(f"delete_library: unexpected reply: {res!r}")


class ObjHistoryClient(ObjDatasetContainerClient[wrappers.History, wrappers.HistoryPreview]):
    """
    Interacts with Galaxy histories.
    """

    CONTAINER_TYPE = wrappers.History
    CONTAINER_PREVIEW_TYPE = wrappers.HistoryPreview

    def create(self, name: Optional[str] = None) -> wrappers.History:
        """
        Create a new Galaxy history, optionally setting its name.

        :rtype: :class:`~.wrappers.History`
        :return: the history just created
        """
        res = self.gi.histories.create_history(name=name)
        hist_info = self._get_dict("create_history", res)
        return self.get(hist_info["id"])

    def list(self, name: Optional[str] = None, deleted: bool = False) -> list[wrappers.History]:
        """
        Get histories owned by the user of this Galaxy instance.

        :type name: str
        :param name: return only histories with this name
        :type deleted: bool
        :param deleted: if ``True``, return histories that have been deleted

        :rtype: list of :class:`~.wrappers.History`
        """
        dicts = self.gi.histories.get_histories(name=name, deleted=deleted)
        return [self.get(_["id"]) for _ in dicts]

    def delete(self, id_: Optional[str] = None, name: Optional[str] = None, purge: bool = False) -> None:
        """
        Delete the history with the given id or name.

        Fails if multiple histories have the same name.

        :type purge: bool
        :param purge: if ``True``, also purge (permanently delete) the history

        .. note::
          The ``purge`` option works only if the Galaxy instance has the
          ``allow_user_dataset_purge`` option set to ``true`` in the
          ``config/galaxy.yml`` configuration file.
        """
        id_ = self._select_id(id_=id_, name=name)
        res = self.gi.histories.delete_history(id_, purge=purge)
        if not isinstance(res, dict):
            raise RuntimeError(f"delete_history: unexpected reply: {res!r}")


class ObjWorkflowClient(ObjClient):
    """
    Interacts with Galaxy workflows.
    """

    def import_new(self, src: Union[str, dict[str, Any]], publish: bool = False) -> wrappers.Workflow:
        """
        Imports a new workflow into Galaxy.

        :type src: dict or str
        :param src: deserialized (dictionary) or serialized (str) JSON
          dump of the workflow (this is normally obtained by exporting
          a workflow from Galaxy).

        :type publish: bool
        :param publish:  if ``True`` the uploaded workflow will be published;
                         otherwise it will be visible only by the user which uploads it (default).

        :rtype: :class:`~.wrappers.Workflow`
        :return: the workflow just imported
        """
        if isinstance(src, dict):
            wf_dict = src
        else:
            try:
                wf_dict = json.loads(src)
            except (TypeError, ValueError):
                raise ValueError(f"src not supported: {src!r}")
        wf_info = self.gi.workflows.import_workflow_dict(wf_dict, publish)
        return self.get(wf_info["id"])

    def import_shared(self, id_: str) -> wrappers.Workflow:
        """
        Imports a shared workflow to the user's space.

        :type id_: str
        :param id_: workflow id

        :rtype: :class:`~.wrappers.Workflow`
        :return: the workflow just imported
        """
        wf_info = self.gi.workflows.import_shared_workflow(id_)
        return self.get(wf_info["id"])

    def get(self, id_: str) -> wrappers.Workflow:
        """
        Retrieve the workflow corresponding to the given id.

        :rtype: :class:`~.wrappers.Workflow`
        :return: the workflow corresponding to ``id_``
        """
        res = self.gi.workflows.show_workflow(id_)
        wf_dict = self._get_dict("show_workflow", res)
        return wrappers.Workflow(wf_dict, gi=self.obj_gi)

    # the 'deleted' option is not available for workflows
    def get_previews(
        self, name: Optional[str] = None, published: bool = False, **kwargs: Any
    ) -> list[wrappers.WorkflowPreview]:
        dicts = self.gi.workflows.get_workflows(name=name, published=published, **kwargs)
        return [wrappers.WorkflowPreview(_, gi=self.obj_gi) for _ in dicts]

    # the 'deleted' option is not available for workflows
    def list(self, name: Optional[str] = None, published: bool = False) -> list[wrappers.Workflow]:
        """
        Get workflows owned by the user of this Galaxy instance.

        :type name: str
        :param name: return only workflows with this name
        :type published: bool
        :param published: if ``True``, return also published workflows

        :rtype: list of :class:`~.wrappers.Workflow`
        """
        dicts = self.gi.workflows.get_workflows(name=name, published=published)
        return [self.get(_["id"]) for _ in dicts]

    def delete(self, id_: Optional[str] = None, name: Optional[str] = None) -> None:
        """
        Delete the workflow with the given id or name.

        Fails if multiple workflows have the specified name.

        .. warning::
          Deleting a workflow is irreversible - all of the data from
          the workflow will be permanently deleted.
        """
        id_ = self._select_id(id_=id_, name=name)
        self.gi.workflows.delete_workflow(id_)


class ObjInvocationClient(ObjClient):
    """
    Interacts with Galaxy Invocations.
    """

    def get(self, id_: str) -> wrappers.Invocation:
        """
        Get an invocation by ID.

        :rtype: Invocation
        :return: invocation object
        """
        inv_dict = self.gi.invocations.show_invocation(id_)
        return wrappers.Invocation(inv_dict, self.obj_gi)

    def get_previews(self, **kwargs: Any) -> list[wrappers.InvocationPreview]:
        """
        Get previews of all invocations.

        :rtype: list of InvocationPreview
        :return: previews of invocations
        """
        inv_list = self.gi.invocations.get_invocations(**kwargs)
        return [wrappers.InvocationPreview(inv_dict, gi=self.obj_gi) for inv_dict in inv_list]

    def list(
        self,
        workflow: Optional[wrappers.Workflow] = None,
        history: Optional[wrappers.History] = None,
        include_terminal: bool = True,
        limit: Optional[int] = None,
        offset: Optional[int] = None,
    ) -> list[wrappers.Invocation]:
        """
        Get full listing of workflow invocations, or select a subset
        by specifying optional arguments for filtering (e.g. a workflow).

        :type workflow: wrappers.Workflow
        :param workflow: Include only invocations associated with
          this workflow

        :type history: wrappers.History
        :param history: Include only invocations associated with
          this history

        :param include_terminal: bool
        :param: Whether to include invocations in terminal state.

        :type limit: int
        :param limit: Maximum number of invocations to return.

        :type offset: int
        :param offset: Number of invocations to skip. Return invocations
           starting from item offset+1.

        :rtype: list of Invocation
        :return: invocation objects
        """
        inv_dict_list = self.gi.invocations.get_invocations(
            workflow_id=workflow.id if workflow else None,
            history_id=history.id if history else None,
            include_terminal=include_terminal,
            limit=limit,
            offset=offset,
            view="element",
            step_details=True,
        )
        return [wrappers.Invocation(inv_dict, self.obj_gi) for inv_dict in inv_dict_list]


class ObjToolClient(ObjClient):
    """
    Interacts with Galaxy tools.
    """

    def get(self, id_: str, io_details: bool = False, link_details: bool = False) -> wrappers.Tool:
        """
        Retrieve the tool corresponding to the given id.

        :type io_details: bool
        :param io_details: if True, get also input and output details

        :type link_details: bool
        :param link_details: if True, get also link details

        :rtype: :class:`~.wrappers.Tool`
        :return: the tool corresponding to ``id_``
        """
        res = self.gi.tools.show_tool(id_, io_details=io_details, link_details=link_details)
        tool_dict = self._get_dict("show_tool", res)
        return wrappers.Tool(tool_dict, gi=self.obj_gi)

    def get_previews(self, name: Optional[str] = None, trackster: bool = False, **kwargs: Any) -> list[wrappers.Tool]:
        """
        Get the list of tools installed on the Galaxy instance.

        :type name: str
        :param name: return only tools with this name

        :type trackster: bool
        :param trackster: if True, only tools that are compatible with
          Trackster are returned

        :rtype: list of :class:`~.wrappers.Tool`
        """
        dicts = self.gi.tools.get_tools(name=name, trackster=trackster, **kwargs)
        return [wrappers.Tool(_, gi=self.obj_gi) for _ in dicts]

    # the 'deleted' option is not available for tools
    def list(self, name: Optional[str] = None, trackster: bool = False) -> list[wrappers.Tool]:
        """
        Get the list of tools installed on the Galaxy instance.

        :type name: str
        :param name: return only tools with this name

        :type trackster: bool
        :param trackster: if True, only tools that are compatible with
          Trackster are returned

        :rtype: list of :class:`~.wrappers.Tool`
        """
        # dicts = self.gi.tools.get_tools(name=name, trackster=trackster)
        # return [self.get(_['id']) for _ in dicts]
        # As of 2015/04/15, GET /api/tools returns also data manager tools for
        # non-admin users, see
        # https://trello.com/c/jyl0cvFP/2633-api-tool-list-filtering-doesn-t-filter-data-managers-for-non-admins
        # Trying to get() a data manager tool would then return a 404 Not Found
        # error.
        # Moreover, the dicts returned by gi.tools.get_tools() are richer than
        # those returned by get(), so make this an alias for get_previews().
        return self.get_previews(name, trackster)


class ObjJobClient(ObjClient):
    """
    Interacts with Galaxy jobs.
    """

    def get(self, id_: str, full_details: bool = False) -> wrappers.Job:
        """
        Retrieve the job corresponding to the given id.

        :type full_details: bool
        :param full_details: if ``True``, return the complete list of details
          for the given job.

        :rtype: :class:`~.wrappers.Job`
        :return: the job corresponding to ``id_``
        """
        res = self.gi.jobs.show_job(id_, full_details)
        job_dict = self._get_dict("show_job", res)
        return wrappers.Job(job_dict, gi=self.obj_gi)

    def get_previews(self, **kwargs: Any) -> list[wrappers.JobPreview]:
        dicts = self.gi.jobs.get_jobs(**kwargs)
        return [wrappers.JobPreview(_, gi=self.obj_gi) for _ in dicts]

    def list(self) -> list[wrappers.Job]:
        """
        Get the list of jobs of the current user.

        :rtype: list of :class:`~.wrappers.Job`
        """
        dicts = self.gi.jobs.get_jobs()
        return [self.get(_["id"]) for _ in dicts]


class ObjDatasetClient(ObjClient):
    """
    Interacts with Galaxy datasets.
    """

    @overload
    def get(self, id_: str, hda_ldda: Literal["hda"] = "hda") -> wrappers.HistoryDatasetAssociation: ...

    @overload
    def get(self, id_: str, hda_ldda: Literal["ldda"]) -> wrappers.LibraryDatasetDatasetAssociation: ...

    def get(self, id_: str, hda_ldda: HdaLdda = "hda") -> wrappers.Dataset:
        """
        Retrieve the dataset corresponding to the given id.

        :type hda_ldda: str
        :param hda_ldda: Whether to show a history dataset ('hda' - the default)
          or library dataset ('ldda')

        :rtype: :class:`~.wrappers.HistoryDatasetAssociation` or :class:`~.wrappers.LibraryDatasetDatasetAssociation`
        :return: the history or library dataset corresponding to ``id_``
        """
        res = self.gi.datasets.show_dataset(id_, hda_ldda=hda_ldda)
        ds_dict = self._get_dict("show_dataset", res)
        if hda_ldda == "hda":
            hist = self.obj_gi.histories.get(ds_dict["history_id"])
            return wrappers.HistoryDatasetAssociation(ds_dict, hist, gi=self.obj_gi)
        elif hda_ldda == "ldda":
            lib = self.obj_gi.libraries.get(ds_dict["parent_library_id"])
            return wrappers.LibraryDatasetDatasetAssociation(ds_dict, lib, gi=self.obj_gi)
        else:
            raise ValueError(f"Unsupported value for hda_ldda: {hda_ldda}")

    def get_previews(self, **kwargs: Any) -> list:
        raise NotImplementedError()

    def list(self) -> list:
        raise NotImplementedError()


class ObjDatasetCollectionClient(ObjClient):
    """
    Interacts with Galaxy dataset collections.
    """

    def get(self, id_: str) -> wrappers.HistoryDatasetCollectionAssociation:
        """
        Retrieve the dataset collection corresponding to the given id.

        :rtype: :class:`~.wrappers.HistoryDatasetCollectionAssociation`
        :return: the history dataset collection corresponding to ``id_``
        """
        res = self.gi.dataset_collections.show_dataset_collection(id_)
        ds_dict = self._get_dict("show_dataset_collection", res)
        hist = self.obj_gi.histories.get(ds_dict["history_id"])
        return wrappers.HistoryDatasetCollectionAssociation(ds_dict, hist, gi=self.obj_gi)

    def get_previews(self, **kwargs: Any) -> list:
        raise NotImplementedError()

    def list(self) -> list:
        raise NotImplementedError()
