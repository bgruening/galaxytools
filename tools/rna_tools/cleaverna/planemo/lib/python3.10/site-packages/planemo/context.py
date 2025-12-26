"""Define the context around Planemo computation.

Abstractions for cross cutting concerns (logging, workspace management,
etc.).
"""

import abc
import logging.config
import os
import shutil
import sys
import traceback
from typing import (
    Any,
    Dict,
    Optional,
    TYPE_CHECKING,
)
from urllib.request import urlopen

from planemo.config import read_global_config

if TYPE_CHECKING:
    from planemo.config import OptionSource


class PlanemoContextInterface(metaclass=abc.ABCMeta):
    """Interface Planemo operations use to access workspace context."""

    @abc.abstractmethod
    def set_option_source(self, param_name, option_source, force=False):
        """Specify how an option was set."""

    @abc.abstractmethod
    def get_option_source(self, param_name):
        """Return OptionSource value indicating how the option was set."""

    @abc.abstractproperty
    def global_config(self):
        """Read Planemo's global configuration."""

    @abc.abstractmethod
    def log(self, msg, *args):
        """Log a message."""

    @abc.abstractmethod
    def vlog(self, msg, *args, **kwds):
        """Log a message only if verbose is enabled."""

    @abc.abstractproperty
    def workspace(self):
        """Create and return Planemo's workspace."""

    @abc.abstractproperty
    def galaxy_profiles_directory(self):
        """Create a return a directory for storing Galaxy profiles."""

    @abc.abstractmethod
    def cache_download(self, url, destination):
        """Use workspace to cache download of ``url``."""


class PlanemoContext(PlanemoContextInterface):
    """Implementation of ``PlanemoContextInterface``"""

    planemo_directory: Optional[str]

    def __init__(self) -> None:
        """Construct a Context object using execution environment."""
        self.home = os.getcwd()
        self._global_config: Optional[Dict] = None
        # Will be set by planemo CLI driver
        self.verbose = False
        self.planemo_config: Optional[str] = None
        self.planemo_directory = None
        self.option_source: Dict[str, "OptionSource"] = {}

    def set_option_source(self, param_name: str, option_source: "OptionSource", force: bool = False) -> None:
        """Specify how an option was set."""
        if not force:
            assert param_name not in self.option_source, f"Option source for [{param_name}] already set"
        self.option_source[param_name] = option_source

    def get_option_source(self, param_name: str) -> "OptionSource":
        """Return OptionSource value indicating how the option was set."""
        assert param_name in self.option_source, f"No option source for [{param_name}]"
        return self.option_source[param_name]

    @property
    def global_config(self) -> Dict[str, Any]:
        """Read Planemo's global configuration.

        As defined most simply by ~/.planemo.yml.
        """
        if self._global_config is None:
            self._global_config = read_global_config(self.planemo_config)
        return self._global_config or {}

    def log(self, msg: str, *args) -> None:
        """Log a message."""
        if args:
            msg %= args
        self._log_message(msg)

    def vlog(self, msg: str, *args, **kwds) -> None:
        """Log a message only if verbose is enabled."""
        if self.verbose:
            self.log(msg, *args)
            if kwds.get("exception", False):
                traceback.print_exc(file=sys.stderr)

    @property
    def workspace(self) -> str:
        """Create and return Planemo's workspace.

        By default this will be ``~/.planemo``.
        """
        if not self.planemo_directory:
            raise Exception("No planemo workspace defined.")
        workspace = self.planemo_directory
        return self._ensure_directory(workspace, "workspace")

    @property
    def galaxy_profiles_directory(self) -> str:
        """Create a return a directory for storing Galaxy profiles."""
        path = os.path.join(self.workspace, "profiles")
        return self._ensure_directory(path, "Galaxy profiles")

    def cache_download(self, url, destination):
        """Use workspace to cache download of ``url``."""
        cache = os.path.join(self.workspace, "cache")
        if not os.path.exists(cache):
            os.makedirs(cache)
        filename = os.path.basename(url)
        cache_destination = os.path.join(cache, filename)
        if not os.path.exists(cache_destination):
            with urlopen(url) as fh:
                content = fh.read()
            if len(content) == 0:
                raise Exception(f"Failed to download [{url}].")
            with open(cache_destination, "wb") as f:
                f.write(content)

        shutil.copy(cache_destination, destination)

    def _ensure_directory(self, path: str, name: str) -> str:
        if not os.path.exists(path):
            os.makedirs(path)
        if not os.path.isdir(path):
            raise Exception(f"Planemo {name} directory [{path}] unavailable.")
        return path

    def _log_message(self, message):
        """Extension point for overriding I/O."""
        print(message)


def configure_standard_planemo_logging(verbose: bool) -> None:
    """Configure Planemo's default logging rules."""
    logging_config = {
        "version": 1,
        "disable_existing_loggers": False,
        "formatters": {
            "verbose": {"format": "%(name)s %(levelname)s %(asctime)s: %(message)s"},
            "simple": {"format": "%(name)s %(levelname)s: %(message)s"},
        },
        "handlers": {
            "console": {
                "level": "DEBUG",
                "class": "logging.StreamHandler",
                "formatter": "simple" if not verbose else "verbose",
            },
        },
        "loggers": {
            # Suppress CWL is beta warning, for Planemo purposes - it is absolutely not.
            "galaxy.tools.parser.factory": {
                "handlers": ["console"],
                "propagate": False,
                "level": "ERROR" if not verbose else "DEBUG",
            },
            "galaxy.tools.deps.commands": {
                "handlers": ["console"],
                "propagate": False,
                "level": "ERROR" if not verbose else "DEBUG",
            },
            "galaxy": {
                "handlers": ["console"],
                "propagate": False,
                "level": "INFO" if not verbose else "DEBUG",
            },
            # @jmchilton
            # I'm fixing up Planemo's lint functionality for CWL and I keep seeing this for the
            # schema metadata stuff (e.g. in the workflows repo). "rdflib.term WARNING:
            # http://schema.org/docs/!DOCTYPE html does not look like a valid URI, trying to
            # serialize this will break.". I'm going to suppress this warning I think, or are the
            # examples wrong and should declare their namespaces differently in some way?
            # @mr-c
            # That particular warning is worth suppressing. A PR to silence it permanently would be very welcome!
            # https://github.com/RDFLib/rdflib/blob/main/rdflib/term.py
            "rdflib.term": {
                "handlers": ["console"],
                "propagate": False,
                "level": "ERROR" if not verbose else "DEBUG",
            },
        },
        "root": {
            "handlers": ["console"],
            "propagate": False,
            "level": "WARNING" if not verbose else "DEBUG",
        },
    }
    logging.config.dictConfig(logging_config)


__all__ = (
    "configure_standard_planemo_logging",
    "PlanemoContextInterface",
    "PlanemoContext",
)
