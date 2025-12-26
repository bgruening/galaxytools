"""Module contianing the :class:`ToilEngine` implementation of :class:`Engine`."""

from typing import (
    Callable,
    List,
    Optional,
)

from planemo import cwl
from planemo.runnable import RunnableType
from .interface import BaseEngine


class ToilEngine(BaseEngine):
    """An :class:`Engine` implementation backed by Toil.

    More information on toil can be found at https://github.com/BD2KGenomics/toil.
    """

    handled_runnable_types = [RunnableType.cwl_tool, RunnableType.cwl_workflow]

    def _run(self, runnables, job_paths, output_collectors: Optional[List[Callable]] = None):
        """Run CWL job using Toil."""
        results = []
        for runnable, job_path in zip(runnables, job_paths):
            results.append(cwl.run_toil(self._ctx, runnable, job_path, **self._kwds))
        return results


__all__ = ("ToilEngine",)
