import json
import os
import sys
from collections.abc import Iterable, Mapping, MutableSequence
from io import BufferedWriter
from typing import IO, TYPE_CHECKING, Any, Callable, Final, Optional, TypeVar, Union

import requests
from rdflib.graph import Graph
from ruamel.yaml.comments import CommentedMap, CommentedSeq
from ruamel.yaml.constructor import RoundTripConstructor
from ruamel.yaml.main import YAML

if TYPE_CHECKING:
    from .fetcher import Fetcher

if sys.version_info >= (3, 11):
    from importlib.resources.abc import Traversable
else:
    from importlib.abc import Traversable

__all__: Final = ["Traversable"]

ContextType = dict[str, Union[dict[str, Any], str, Iterable[str]]]
DocumentType = TypeVar("DocumentType", CommentedSeq, CommentedMap)
DocumentOrStrType = TypeVar("DocumentOrStrType", CommentedSeq, CommentedMap, str)
FieldType = TypeVar("FieldType", str, CommentedSeq, CommentedMap)
MandatoryResolveType = Union[int, float, str, CommentedMap, CommentedSeq]
ResolveType = Optional[MandatoryResolveType]
ResolvedRefType = tuple[ResolveType, CommentedMap]
IdxResultType = Union[CommentedMap, CommentedSeq, str, None]
IdxType = dict[str, IdxResultType]
CacheType = dict[str, Union[str, Graph, bool]]
FetcherCallableType = Callable[[CacheType, requests.sessions.Session], "Fetcher"]
AttachmentsType = Callable[[Union[CommentedMap, CommentedSeq]], bool]


def add_dictlist(di: dict[Any, Any], key: Any, val: Any) -> None:
    """Manage element insertion in dicts of lists."""
    if key not in di:
        di[key] = []
    di[key].append(val)


def aslist(thing: Any) -> MutableSequence[Any]:
    """
    Wrap single items and lists.

    Return lists unchanged.
    """
    if isinstance(thing, MutableSequence):
        return thing
    return [thing]


def flatten(thing: Any, ltypes: Any = (list, tuple)) -> Any:
    """Flatten lists recursively."""
    # http://rightfootin.blogspot.com/2006/09/more-on-python-flatten.html
    if thing is None:
        return []
    if not isinstance(thing, ltypes):
        return [thing]

    ltype: Final = type(thing)
    lst: Final = list(thing)
    i = 0
    while i < len(lst):
        while isinstance(lst[i], ltypes):
            if not lst[i]:
                lst.pop(i)
                i -= 1
                break
            lst[i : i + 1] = lst[i]
        i += 1
    return ltype(lst)


def onWindows() -> bool:
    """Check if Python is running on Windows OS."""
    return os.name == "nt"


def convert_to_dict(j4: Any) -> Any:
    """Convert generic Mapping objects to dicts recursively."""
    if isinstance(j4, Mapping):
        return {k: convert_to_dict(v) for k, v in j4.items()}
    if isinstance(j4, MutableSequence):
        return [convert_to_dict(v) for v in j4]
    return j4


def json_dump(obj: Any, fp: IO[str], **kwargs: Any) -> None:
    """Force use of unicode."""
    json.dump(convert_to_dict(obj), fp, **kwargs)


def json_dumps(
    obj: Any,
    **kwargs: Any,
) -> str:
    """Force use of unicode."""
    return json.dumps(convert_to_dict(obj), **kwargs)


def stdout() -> BufferedWriter:
    """Build a replacement for sys.stdout that allow for writing binary data."""
    return os.fdopen(sys.stdout.fileno(), "wb", closefd=False)


class _RoundTripNoTimeStampConstructor(RoundTripConstructor):
    def construct_yaml_timestamp(self: Any, node: Any, values: Any = None) -> Any:
        return node.value


_RoundTripNoTimeStampConstructor.add_constructor(
    "tag:yaml.org,2002:timestamp",
    _RoundTripNoTimeStampConstructor.construct_yaml_timestamp,
)

# mypy: no-warn-unused-ignores


def yaml_no_ts() -> YAML:
    """
    Get a YAML loader that won't parse timestamps into datetime objects.

    Such datetime objects can't be easily dumped into JSON.
    """
    yaml: Final = YAML(typ="rt")
    yaml.preserve_quotes = True  # type: ignore
    yaml.Constructor = _RoundTripNoTimeStampConstructor
    return yaml
