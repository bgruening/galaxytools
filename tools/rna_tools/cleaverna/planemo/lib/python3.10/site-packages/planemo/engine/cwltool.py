"""Module contianing the :class:`CwlToolEngine` implementation of :class:`Engine`."""

from typing import (
    Callable,
    List,
    Optional,
)

from planemo import cwl
from planemo.runnable import RunnableType
from .interface import BaseEngine


class CwlToolEngine(BaseEngine):
    """An :class:`Engine` implementation backed by cwltool.

    More information on cwltool can be found at https://github.com/common-workflow-language/cwltool.
    """

    handled_runnable_types = [RunnableType.cwl_tool, RunnableType.cwl_workflow]

    def _run(self, runnables, job_paths, output_collectors: Optional[List[Callable]] = None):
        """Run CWL job using cwltool."""
        results = []
        for runnable, job_path in zip(runnables, job_paths):
            results.append(cwl.run_cwltool(self._ctx, runnable, job_path, **self._kwds))
        return results


__all__ = ("CwlToolEngine",)
