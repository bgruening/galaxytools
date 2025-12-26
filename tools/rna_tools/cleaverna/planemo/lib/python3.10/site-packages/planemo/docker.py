"""Docker utilities for planemo.

Built on Galaxy abstractions in :mod:`galaxy.tools.deps.dockerfiles` and
:mod:`galaxy.tools.deps.docker_util`.
"""

from galaxy.tool_util.deps.dockerfiles import docker_host_args

__all__ = ("docker_host_args",)
