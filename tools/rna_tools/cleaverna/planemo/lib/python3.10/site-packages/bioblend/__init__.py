import contextlib
import logging
import logging.config
import os
import time
from typing import (
    Callable,
    Optional,
    TypeVar,
    Union,
)

from bioblend.config import (
    BioBlendConfigLocations,
    Config,
)

# Current version of the library
__version__ = "1.7.0"

# default chunk size (in bytes) for reading remote data
try:
    import resource

    CHUNK_SIZE = resource.getpagesize()
except Exception:
    CHUNK_SIZE = 4096


config = Config()


def get_version() -> str:
    """
    Returns a string with the current version of the library (e.g., "0.2.0")
    """
    return __version__


def init_logging() -> None:
    """
    Initialize BioBlend's logging from a configuration file.
    """
    for config_file in BioBlendConfigLocations:
        with contextlib.suppress(Exception):
            logging.config.fileConfig(os.path.expanduser(config_file))


class NullHandler(logging.Handler):
    def emit(self, record: logging.LogRecord) -> None:
        pass


# By default, do not force any logging by the library. If you want to see the
# log messages in your scripts, add the following to the top of your script:
#   import logging
#   logging.basicConfig(filename="bioblend.log", level=logging.DEBUG)
default_format_string = "%(asctime)s %(name)s [%(levelname)s]: %(message)s"
log = logging.getLogger("bioblend")
log.addHandler(NullHandler())
init_logging()

# Convenience functions to set logging to a particular file or stream
# To enable either of these, simply add the following at the top of a
# bioblend module:
#   import bioblend
#   bioblend.set_stream_logger(__name__)


def set_file_logger(
    name: str, filepath: str, level: Union[int, str] = logging.INFO, format_string: Optional[str] = None
) -> None:
    global log
    if not format_string:
        format_string = default_format_string
    logger = logging.getLogger(name)
    logger.setLevel(level)
    fh = logging.FileHandler(filepath)
    fh.setLevel(level)
    formatter = logging.Formatter(format_string)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    log = logger


def set_stream_logger(name: str, level: Union[int, str] = logging.DEBUG, format_string: Optional[str] = None) -> None:
    global log
    if not format_string:
        format_string = default_format_string
    logger = logging.getLogger(name)
    logger.setLevel(level)
    fh = logging.StreamHandler()
    fh.setLevel(level)
    formatter = logging.Formatter(format_string)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    log = logger


class ConnectionError(Exception):
    """
    An exception class that is raised when unexpected HTTP responses come back.

    Should make it easier to debug when strange HTTP things happen such as a
    proxy server getting in the way of the request etc.
    @see: body attribute to see the content of the http response
    """

    def __init__(  # noqa: B042  # https://github.com/PyCQA/flake8-bugbear/issues/525
        self, message: str, body: Optional[Union[bytes, str]] = None, status_code: Optional[int] = None
    ) -> None:
        super().__init__(message)
        self.body = body
        self.status_code = status_code

    def __str__(self) -> str:
        return f"{self.args[0]}: {self.body!s}"


class TimeoutException(Exception):
    pass


class NotReady(Exception):
    pass


T = TypeVar("T")


def wait_on(func: Callable[[], T], maxwait: float = 60, interval: float = 3) -> T:
    """
    Wait until a function returns without raising a NotReady exception

    :param func: function to wait on. It should accept no parameters.

    :param maxwait: Total time (in seconds) to wait for the function to return
      without raising a NotReady exception. After this time, a
      ``TimeoutException`` will be raised.

    :param interval: Time (in seconds) to wait between 2 consecutive checks.
    """
    assert maxwait >= 0
    assert interval > 0

    time_left = maxwait
    while True:
        try:
            return func()
        except NotReady as e:
            if time_left > 0:
                log.info("%s. Will wait %s more s", e, time_left)
                time.sleep(min(time_left, interval))
                time_left -= interval
            else:
                raise TimeoutException(f"{e} after {maxwait} s")
