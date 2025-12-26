import logging
from typing import (
    Any,
    Optional,
    TYPE_CHECKING,
    Union,
)

from bioblend import (
    CHUNK_SIZE,
    NotReady,
    TimeoutException,
    wait_on,
)
from bioblend.galaxy.client import Client
from bioblend.galaxy.datasets import TERMINAL_STATES

if TYPE_CHECKING:
    from bioblend.galaxy import GalaxyInstance

log = logging.getLogger(__name__)


class HasElements:
    def __init__(
        self,
        name: str,
        type: str = "list",
        elements: Optional[Union[list[Union["CollectionElement", "SimpleElement"]], dict[str, Any]]] = None,
    ) -> None:
        self.name = name
        self.type = type
        if isinstance(elements, dict):
            self.elements: list[Union[CollectionElement, SimpleElement]] = [
                HistoryDatasetElement(name=key, id=value) for key, value in elements.values()
            ]
        elif elements:
            self.elements = elements

    def add(self, element: Union["CollectionElement", "SimpleElement"]) -> "HasElements":
        self.elements.append(element)
        return self


class CollectionDescription(HasElements):
    def to_dict(self) -> dict[str, Union[str, list]]:
        return {
            "name": self.name,
            "collection_type": self.type,
            "element_identifiers": [e.to_dict() for e in self.elements],
        }


class CollectionElement(HasElements):
    def to_dict(self) -> dict[str, Union[str, list]]:
        return {
            "src": "new_collection",
            "name": self.name,
            "collection_type": self.type,
            "element_identifiers": [e.to_dict() for e in self.elements],
        }


class SimpleElement:
    def __init__(self, value: dict[str, str]) -> None:
        self.value = value

    def to_dict(self) -> dict[str, str]:
        return self.value


class HistoryDatasetElement(SimpleElement):
    def __init__(self, name: str, id: str) -> None:
        super().__init__(
            {
                "name": name,
                "src": "hda",
                "id": id,
            }
        )


class HistoryDatasetCollectionElement(SimpleElement):
    def __init__(self, name: str, id: str) -> None:
        super().__init__(
            {
                "name": name,
                "src": "hdca",
                "id": id,
            }
        )


class LibraryDatasetElement(SimpleElement):
    def __init__(self, name: str, id: str) -> None:
        super().__init__(
            {
                "name": name,
                "src": "ldda",
                "id": id,
            }
        )


class DatasetCollectionClient(Client):
    gi: "GalaxyInstance"
    module = "dataset_collections"

    def __init__(self, galaxy_instance: "GalaxyInstance") -> None:
        super().__init__(galaxy_instance)

    def show_dataset_collection(self, dataset_collection_id: str, instance_type: str = "history") -> dict[str, Any]:
        """
        Get details of a given dataset collection of the current user

        :type dataset_collection_id: str
        :param dataset_collection_id: dataset collection ID

        :type instance_type: str
        :param instance_type: instance type of the collection - 'history' or 'library'

        :rtype: dict
        :return: element view of the dataset collection
        """
        params = {
            "instance_type": instance_type,
        }
        url = self._make_url(module_id=dataset_collection_id)
        return self._get(id=dataset_collection_id, url=url, params=params)

    def download_dataset_collection(self, dataset_collection_id: str, file_path: str) -> dict[str, Any]:
        """
        Download a history dataset collection as an archive.

        :type dataset_collection_id: str
        :param dataset_collection_id: Encoded dataset collection ID

        :type file_path: str
        :param file_path: The path to which the archive will be downloaded

        :rtype: dict
        :return: Information about the downloaded archive.

        .. note::
          This method downloads a ``zip`` archive for Galaxy 21.01 and later.
          For earlier versions of Galaxy this method downloads a ``tgz`` archive.
        """
        url = self._make_url(module_id=dataset_collection_id) + "/download"
        r = self.gi.make_get_request(url, stream=True)
        r.raise_for_status()

        archive_type = "zip" if self.gi.config.get_version()["version_major"] >= "21.01" else "tgz"

        with open(file_path, "wb") as fp:
            for chunk in r.iter_content(chunk_size=CHUNK_SIZE):
                if chunk:
                    fp.write(chunk)

        return {"file_path": file_path, "archive_type": archive_type}

    def wait_for_dataset_collection(
        self,
        dataset_collection_id: str,
        maxwait: float = 12000,
        interval: float = 3,
        proportion_complete: float = 1.0,
        check: bool = True,
    ) -> dict[str, Any]:
        """
        Wait until all or a specified proportion of elements of a dataset
        collection are in a terminal state.

        :type dataset_collection_id: str
        :param dataset_collection_id: dataset collection ID

        :type maxwait: float
        :param maxwait: Total time (in seconds) to wait for the dataset
          states in the dataset collection to become terminal. After this time,
          a ``TimeoutException`` will be raised.

        :type interval: float
        :param interval: Time (in seconds) to wait between two consecutive checks.

        :type proportion_complete: float
        :param proportion_complete: Proportion of elements in this collection
          that have to be in a terminal state for this method to return. Must
          be a number between 0 and 1. For example: if the dataset collection
          contains 2 elements, and proportion_complete=0.5 is specified, then
          wait_for_dataset_collection will return as soon as 1 of the 2
          datasets is in a terminal state. Default is 1, i.e. all elements must
          complete.

        :type check: bool
        :param check: Whether to check if all the terminal states of datasets
          in the dataset collection are 'ok'. This will raise an Exception if
          a dataset is in a terminal state other than 'ok'.

        :rtype: dict
        :return: Details of the given dataset collection.
        """
        assert 0 <= proportion_complete <= 1

        def check_and_get_dataset_collection() -> dict[str, Any]:
            dataset_collection = self.show_dataset_collection(dataset_collection_id)
            states = [elem["object"]["state"] for elem in dataset_collection["elements"]]
            terminal_states = [state for state in states if state in TERMINAL_STATES]
            if set(terminal_states) not in [{"ok"}, set()]:
                raise Exception(
                    f"Dataset collection {dataset_collection_id} contains elements in the "
                    f"following non-ok terminal states: {', '.join(set(terminal_states) - {'ok'})}"
                )
            proportion = len(terminal_states) / len(states)
            if proportion >= proportion_complete:
                return dataset_collection
            raise NotReady(
                f"The dataset collection {dataset_collection_id} has only {proportion * 100}% of datasets in a terminal state"
            )

        return wait_on(check_and_get_dataset_collection, maxwait=maxwait, interval=interval)


# Unused, for backward compatibility
DatasetCollectionTimeoutException = TimeoutException


__all__ = (
    "CollectionDescription",
    "CollectionElement",
    "DatasetCollectionClient",
    "HistoryDatasetElement",
    "HistoryDatasetCollectionElement",
    "LibraryDatasetElement",
)
