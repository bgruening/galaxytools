"""
A representation of a Galaxy instance based on oo wrappers.
"""

import time
from collections.abc import Iterable
from typing import Optional

import bioblend
import bioblend.galaxy
from bioblend.galaxy.datasets import TERMINAL_STATES
from . import (
    client,
    wrappers,
)


def _get_error_info(dataset: wrappers.Dataset) -> str:
    msg = dataset.id
    try:
        msg += f" ({dataset.name}): {dataset.misc_info}"
    except Exception:  # avoid 'error while generating an error report'
        msg += ": error"
    return msg


class GalaxyInstance:
    """
    A representation of an instance of Galaxy, identified by a URL and
    a user's API key.

    :type url: str
    :param url: a FQDN or IP for a given instance of Galaxy. For example:
      ``http://127.0.0.1:8080``

    :type api_key: str
    :param api_key: user's API key for the given instance of Galaxy, obtained
      from the Galaxy web UI.

    This is actually a factory class which instantiates the entity-specific
    clients.

    Example: get a list of all histories for a user with API key 'foo'::

      from bioblend.galaxy.objects import GalaxyInstance
      gi = GalaxyInstance('http://127.0.0.1:8080', api_key='foo')
      histories = gi.histories.list()
    """

    def __init__(
        self,
        url: str,
        api_key: Optional[str] = None,
        email: Optional[str] = None,
        password: Optional[str] = None,
        *,
        token: Optional[str] = None,
        verify: bool = True,
        user_agent: Optional[str] = None,
    ) -> None:
        self.gi = bioblend.galaxy.GalaxyInstance(
            url, key=api_key, email=email, password=password, token=token, verify=verify, user_agent=user_agent
        )
        self.log = bioblend.log
        self.datasets = client.ObjDatasetClient(self)
        self.dataset_collections = client.ObjDatasetCollectionClient(self)
        self.histories = client.ObjHistoryClient(self)
        self.libraries = client.ObjLibraryClient(self)
        self.workflows = client.ObjWorkflowClient(self)
        self.invocations = client.ObjInvocationClient(self)
        self.tools = client.ObjToolClient(self)
        self.jobs = client.ObjJobClient(self)

    def _wait_datasets(
        self, datasets: Iterable[wrappers.Dataset], polling_interval: float, break_on_error: bool = True
    ) -> None:
        """
        Wait for datasets to come out of the pending states.

        :type datasets: :class:`~collections.Iterable` of
          :class:`~.wrappers.Dataset`
        :param datasets: datasets

        :type polling_interval: float
        :param polling_interval: polling interval in seconds

        :type break_on_error: bool
        :param break_on_error: if ``True``, raise a RuntimeError exception as
          soon as at least one of the datasets is in the 'error' state.

        .. warning::

          This is a blocking operation that can take a very long time.
          Also, note that this method does not return anything;
          however, each input dataset is refreshed (possibly multiple
          times) during the execution.
        """

        def poll(ds_list: Iterable[wrappers.Dataset]) -> list[wrappers.Dataset]:
            pending = []
            for ds in ds_list:
                ds.refresh()
                if break_on_error and ds.state == "error":
                    raise RuntimeError(_get_error_info(ds))
                if not ds.state:
                    self.log.warning("Dataset %s has an empty state", ds.id)
                elif ds.state not in TERMINAL_STATES:
                    self.log.info("Dataset %s is in non-terminal state %s", ds.id, ds.state)
                    pending.append(ds)
            return pending

        self.log.info("Waiting for datasets")
        while datasets:
            datasets = poll(datasets)
            time.sleep(polling_interval)
