from typing import List

NO_GLOB_2 = "glob2 library unavailabile, please install with pip install glob2."

try:
    from glob2 import glob as _glob
except ImportError:
    _glob = None


def glob(*args, **kwds) -> List[str]:
    if _glob is None:
        raise Exception(NO_GLOB_2)
    return _glob(*args, **kwds)


__all__ = ("glob",)
