"""Template code used by python_codegen.py."""

import copy
import logging
import os
import pathlib
import tempfile
import uuid as _uuid__  # pylint: disable=unused-import # noqa: F401
import xml.sax  # nosec
from abc import ABC, abstractmethod
from collections.abc import MutableMapping, MutableSequence, Sequence
from io import StringIO
from itertools import chain
from typing import Any, Final, Optional, Union, cast
from urllib.parse import quote, urldefrag, urlparse, urlsplit, urlunsplit
from urllib.request import pathname2url

from rdflib import Graph
from rdflib.plugins.parsers.notation3 import BadSyntax
from ruamel.yaml.comments import CommentedMap

from schema_salad.exceptions import SchemaSaladException, ValidationException
from schema_salad.fetcher import DefaultFetcher, Fetcher, MemoryCachingFetcher
from schema_salad.sourceline import SourceLine, add_lc_filename
from schema_salad.utils import CacheType, yaml_no_ts  # requires schema-salad v8.2+

_vocab: dict[str, str] = {}
_rvocab: dict[str, str] = {}

_logger: Final = logging.getLogger("salad")


IdxType = MutableMapping[str, tuple[Any, "LoadingOptions"]]


class LoadingOptions:
    idx: Final[IdxType]
    fileuri: Final[Optional[str]]
    baseuri: Final[str]
    namespaces: Final[MutableMapping[str, str]]
    schemas: Final[MutableSequence[str]]
    original_doc: Final[Optional[Any]]
    addl_metadata: Final[MutableMapping[str, Any]]
    fetcher: Final[Fetcher]
    vocab: Final[dict[str, str]]
    rvocab: Final[dict[str, str]]
    cache: Final[CacheType]
    imports: Final[list[str]]
    includes: Final[list[str]]
    no_link_check: Final[Optional[bool]]
    container: Final[Optional[str]]

    def __init__(
        self,
        fetcher: Optional[Fetcher] = None,
        namespaces: Optional[dict[str, str]] = None,
        schemas: Optional[list[str]] = None,
        fileuri: Optional[str] = None,
        copyfrom: Optional["LoadingOptions"] = None,
        original_doc: Optional[Any] = None,
        addl_metadata: Optional[dict[str, str]] = None,
        baseuri: Optional[str] = None,
        idx: Optional[IdxType] = None,
        imports: Optional[list[str]] = None,
        includes: Optional[list[str]] = None,
        no_link_check: Optional[bool] = None,
        container: Optional[str] = None,
    ) -> None:
        """Create a LoadingOptions object."""
        self.original_doc = original_doc

        if idx is not None:
            temp_idx = idx
        else:
            temp_idx = copyfrom.idx if copyfrom is not None else {}
        self.idx = temp_idx

        if fileuri is not None:
            temp_fileuri: Optional[str] = fileuri
        else:
            temp_fileuri = copyfrom.fileuri if copyfrom is not None else None
        self.fileuri = temp_fileuri

        if baseuri is not None:
            temp_baseuri = baseuri
        else:
            temp_baseuri = copyfrom.baseuri if copyfrom is not None else ""
        self.baseuri = temp_baseuri

        if namespaces is not None:
            temp_namespaces: MutableMapping[str, str] = namespaces
        else:
            temp_namespaces = copyfrom.namespaces if copyfrom is not None else {}
        self.namespaces = temp_namespaces

        if schemas is not None:
            temp_schemas: MutableSequence[str] = schemas
        else:
            temp_schemas = copyfrom.schemas if copyfrom is not None else []
        self.schemas = temp_schemas

        if addl_metadata is not None:
            temp_addl_metadata: MutableMapping[str, Any] = addl_metadata
        else:
            temp_addl_metadata = copyfrom.addl_metadata if copyfrom is not None else {}
        self.addl_metadata = temp_addl_metadata

        if imports is not None:
            temp_imports = imports
        else:
            temp_imports = copyfrom.imports if copyfrom is not None else []
        self.imports = temp_imports

        if includes is not None:
            temp_includes = includes
        else:
            temp_includes = copyfrom.includes if copyfrom is not None else []
        self.includes = temp_includes

        if no_link_check is not None:
            temp_no_link_check: Optional[bool] = no_link_check
        else:
            temp_no_link_check = copyfrom.no_link_check if copyfrom is not None else False
        self.no_link_check = temp_no_link_check

        if container is not None:
            temp_container: Optional[str] = container
        else:
            temp_container = copyfrom.container if copyfrom is not None else None
        self.container = temp_container

        if fetcher is not None:
            temp_fetcher = fetcher
        elif copyfrom is not None:
            temp_fetcher = copyfrom.fetcher
        else:
            import requests
            from cachecontrol.caches import SeparateBodyFileCache
            from cachecontrol.wrapper import CacheControl

            root = pathlib.Path(os.environ.get("HOME", tempfile.gettempdir()))
            session = CacheControl(
                requests.Session(),
                cache=SeparateBodyFileCache(root / ".cache" / "salad"),
            )
            temp_fetcher = DefaultFetcher({}, session)
        self.fetcher = temp_fetcher

        self.cache = self.fetcher.cache if isinstance(self.fetcher, MemoryCachingFetcher) else {}

        if self.namespaces != {}:
            temp_vocab = _vocab.copy()
            temp_rvocab = _rvocab.copy()
            for k, v in self.namespaces.items():
                temp_vocab[k] = v
                temp_rvocab[v] = k
        else:
            temp_vocab = _vocab
            temp_rvocab = _rvocab
        self.vocab = temp_vocab
        self.rvocab = temp_rvocab

    @property
    def graph(self) -> Graph:
        """Generate a merged rdflib.Graph from all entries in self.schemas."""
        graph = Graph()
        if not self.schemas:
            return graph
        key: Final = str(hash(tuple(self.schemas)))
        if key in self.cache:
            return cast(Graph, self.cache[key])
        for schema in self.schemas:
            fetchurl = (
                self.fetcher.urljoin(self.fileuri, schema)
                if self.fileuri is not None
                else pathlib.Path(schema).resolve().as_uri()
            )
            if fetchurl not in self.cache or self.cache[fetchurl] is True:
                _logger.debug("Getting external schema %s", fetchurl)
                try:
                    content = self.fetcher.fetch_text(fetchurl)
                except Exception as e:
                    _logger.warning("Could not load extension schema %s: %s", fetchurl, str(e))
                    continue
                newGraph = Graph()
                err_msg = "unknown error"
                for fmt in ["xml", "turtle"]:
                    try:
                        newGraph.parse(data=content, format=fmt, publicID=str(fetchurl))
                        self.cache[fetchurl] = newGraph
                        graph += newGraph
                        break
                    except (xml.sax.SAXParseException, TypeError, BadSyntax) as e:
                        err_msg = str(e)
                else:
                    _logger.warning("Could not load extension schema %s: %s", fetchurl, err_msg)
        self.cache[key] = graph
        return graph


class Saveable(ABC):
    """Mark classes than have a save() and fromDoc() function."""

    @classmethod
    @abstractmethod
    def fromDoc(
        cls,
        _doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
    ) -> "Saveable":
        """Construct this object from the result of yaml.load()."""

    @abstractmethod
    def save(
        self, top: bool = False, base_url: str = "", relative_uris: bool = True
    ) -> dict[str, Any]:
        """Convert this object to a JSON/YAML friendly dictionary."""


def load_field(
    val: Union[str, dict[str, str]],
    fieldtype: "_Loader",
    baseuri: str,
    loadingOptions: LoadingOptions,
    lc: Optional[list[Any]] = None,
) -> Any:
    """Load field."""
    if isinstance(val, MutableMapping):
        if "$import" in val:
            if loadingOptions.fileuri is None:
                raise SchemaSaladException("Cannot load $import without fileuri")
            url1: Final = loadingOptions.fetcher.urljoin(loadingOptions.fileuri, val["$import"])
            result, metadata = _document_load_by_url(
                fieldtype,
                url1,
                loadingOptions,
            )
            loadingOptions.imports.append(url1)
            return result
        if "$include" in val:
            if loadingOptions.fileuri is None:
                raise SchemaSaladException("Cannot load $import without fileuri")
            url2: Final = loadingOptions.fetcher.urljoin(loadingOptions.fileuri, val["$include"])
            val = loadingOptions.fetcher.fetch_text(url2)
            loadingOptions.includes.append(url2)
    return fieldtype.load(val, baseuri, loadingOptions, lc=lc)


save_type = Optional[Union[MutableMapping[str, Any], MutableSequence[Any], int, float, bool, str]]


def extract_type(val_type: type[Any]) -> str:
    """Take a type of value, and extracts the value as a string."""
    val_str: Final = str(val_type)
    return val_str.split("'")[1]


def convert_typing(val_type: str) -> str:
    """Normalize type names to schema-salad types."""
    if "None" in val_type:
        return "null"
    if "CommentedSeq" in val_type or "list" in val_type:
        return "array"
    if "CommentedMap" in val_type or "dict" in val_type:
        return "object"
    if "False" in val_type or "True" in val_type:
        return "boolean"
    return val_type


def parse_errors(error_message: str) -> tuple[str, str, str]:
    """Parse error messages from several loaders into one error message."""
    if not error_message.startswith("Expected"):
        return error_message, "", ""
    vals: Final = error_message.split("\n")
    if len(vals) == 1:
        return error_message, "", ""
    types1: Final = set()
    for val in vals:
        individual_vals = val.split(" ")
        if val == "":
            continue
        if individual_vals[1] == "one":
            individual_vals = val.split("(")[1].split(",")
            for t in individual_vals:
                types1.add(t.strip(" ").strip(")\n"))
        elif individual_vals[2] == "<class":
            types1.add(individual_vals[3].strip(">").replace("'", ""))
        elif individual_vals[0] == "Value":
            types1.add(individual_vals[-1].strip("."))
        else:
            types1.add(individual_vals[1].replace(",", ""))
    types2: Final = {val for val in types1 if val != "NoneType"}
    if "str" in types2:
        types3 = {convert_typing(val) for val in types2 if "'" not in val}
    else:
        types3 = types2
    to_print = ""
    for val in types3:
        if "'" in val:
            to_print = "value" if len(types3) == 1 else "values"

    if to_print == "":
        to_print = "type" if len(types3) == 1 else "types"

    verb_tensage: Final = "is" if len(types3) == 1 else "are"

    return str(types3).replace("{", "(").replace("}", ")").replace("'", ""), to_print, verb_tensage


def save(
    val: Any,
    top: bool = True,
    base_url: str = "",
    relative_uris: bool = True,
) -> save_type:
    if isinstance(val, Saveable):
        return val.save(top=top, base_url=base_url, relative_uris=relative_uris)
    if isinstance(val, MutableSequence):
        return [save(v, top=False, base_url=base_url, relative_uris=relative_uris) for v in val]
    if isinstance(val, MutableMapping):
        newdict: Final = {}
        for key in val:
            newdict[key] = save(val[key], top=False, base_url=base_url, relative_uris=relative_uris)
        return newdict
    if val is None or isinstance(val, (int, float, bool, str)):
        return val
    raise Exception("Not Saveable: %s" % type(val))


def save_with_metadata(
    val: Any,
    valLoadingOpts: LoadingOptions,
    top: bool = True,
    base_url: str = "",
    relative_uris: bool = True,
) -> save_type:
    """Save and set $namespaces, $schemas, $base and any other metadata fields at the top level."""
    saved_val: Final = save(val, top, base_url, relative_uris)
    newdict: MutableMapping[str, Any] = {}
    if isinstance(saved_val, MutableSequence):
        newdict = {"$graph": saved_val}
    elif isinstance(saved_val, MutableMapping):
        newdict = saved_val

    if valLoadingOpts.namespaces:
        newdict["$namespaces"] = valLoadingOpts.namespaces
    if valLoadingOpts.schemas:
        newdict["$schemas"] = valLoadingOpts.schemas
    if valLoadingOpts.baseuri:
        newdict["$base"] = valLoadingOpts.baseuri
    for k, v in valLoadingOpts.addl_metadata.items():
        if k not in newdict:
            newdict[k] = v

    return newdict


def expand_url(
    url: str,
    base_url: str,
    loadingOptions: LoadingOptions,
    scoped_id: bool = False,
    vocab_term: bool = False,
    scoped_ref: Optional[int] = None,
) -> str:
    if url in ("@id", "@type"):
        return url

    if vocab_term and url in loadingOptions.vocab:
        return url

    if bool(loadingOptions.vocab) and ":" in url:
        prefix: Final = url.split(":")[0]
        if prefix in loadingOptions.vocab:
            url = loadingOptions.vocab[prefix] + url[len(prefix) + 1 :]

    split1: Final = urlsplit(url)

    if (
        (bool(split1.scheme) and split1.scheme in loadingOptions.fetcher.supported_schemes())
        or url.startswith("$(")
        or url.startswith("${")
    ):
        pass
    elif scoped_id and not bool(split1.fragment):
        splitbase1: Final = urlsplit(base_url)
        frg: str
        if bool(splitbase1.fragment):
            frg = splitbase1.fragment + "/" + split1.path
        else:
            frg = split1.path
        pt: Final = splitbase1.path if splitbase1.path != "" else "/"
        url = urlunsplit((splitbase1.scheme, splitbase1.netloc, pt, splitbase1.query, frg))
    elif scoped_ref is not None and not bool(split1.fragment):
        splitbase2: Final = urlsplit(base_url)
        sp = splitbase2.fragment.split("/")
        n = scoped_ref
        while n > 0 and len(sp) > 0:
            sp.pop()
            n -= 1
        sp.append(url)
        url = urlunsplit(
            (
                splitbase2.scheme,
                splitbase2.netloc,
                splitbase2.path,
                splitbase2.query,
                "/".join(sp),
            )
        )
    else:
        url = loadingOptions.fetcher.urljoin(base_url, url)

    if vocab_term:
        split2: Final = urlsplit(url)
        if bool(split2.scheme):
            if url in loadingOptions.rvocab:
                return loadingOptions.rvocab[url]
        else:
            raise ValidationException(f"Term {url!r} not in vocabulary")

    return url


class _Loader:
    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        pass


class _AnyLoader(_Loader):
    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if doc is not None:
            return doc
        raise ValidationException("Expected non-null")


class _PrimitiveLoader(_Loader):
    def __init__(self, tp: Union[type, tuple[type[str], type[str]]]) -> None:
        self.tp: Final = tp

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if not isinstance(doc, self.tp):
            raise ValidationException(f"Expected a {self.tp} but got {doc.__class__.__name__}")
        return doc

    def __repr__(self) -> str:
        return str(self.tp)


class _ArrayLoader(_Loader):
    def __init__(self, items: _Loader) -> None:
        self.items: Final = items

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if not isinstance(doc, MutableSequence):
            raise ValidationException(
                f"Value is a {convert_typing(extract_type(type(doc)))}, "
                f"but valid type for this field is an array."
            )
        r: Final[list[Any]] = []
        errors: Final[list[SchemaSaladException]] = []
        fields: Final[list[str]] = []
        for i in range(0, len(doc)):
            try:
                lf = load_field(
                    doc[i], _UnionLoader([self, self.items]), baseuri, loadingOptions, lc=lc
                )
                flatten = loadingOptions.container != "@list"
                if flatten and isinstance(lf, MutableSequence):
                    r.extend(lf)
                else:
                    r.append(lf)

                if isinstance(doc[i], CommentedMap):
                    if doc[i].get("id") is not None:
                        if doc[i].get("id") in fields:
                            errors.append(
                                ValidationException(
                                    f"Duplicate field {doc[i].get('id')!r}",
                                    SourceLine(doc[i], "id", str),
                                    [],
                                )
                            )
                        else:
                            fields.append(doc[i].get("id"))

            except ValidationException as e:
                e = ValidationException(
                    "array item is invalid because", SourceLine(doc, i, str), [e]
                )
                errors.append(e)
        if errors:
            raise ValidationException("", None, errors)
        return r

    def __repr__(self) -> str:
        return f"array<{self.items}>"


class _MapLoader(_Loader):
    def __init__(
        self,
        values: _Loader,
        name: Optional[str] = None,
        container: Optional[str] = None,
        no_link_check: Optional[bool] = None,
    ) -> None:
        self.values: Final = values
        self.name: Final = name
        self.container: Final = container
        self.no_link_check: Final = no_link_check

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if not isinstance(doc, MutableMapping):
            raise ValidationException(f"Expected a map, was {type(doc)}")
        if self.container is not None or self.no_link_check is not None:
            loadingOptions = LoadingOptions(
                copyfrom=loadingOptions, container=self.container, no_link_check=self.no_link_check
            )
        r: Final[dict[str, Any]] = {}
        errors: Final[list[SchemaSaladException]] = []
        for k, v in doc.items():
            try:
                lf = load_field(v, self.values, baseuri, loadingOptions, lc)
                r[k] = lf
            except ValidationException as e:
                errors.append(e.with_sourceline(SourceLine(doc, k, str)))
        if errors:
            raise ValidationException("", None, errors)
        return r

    def __repr__(self) -> str:
        return self.name if self.name is not None else f"map<string, {self.values}>"


class _EnumLoader(_Loader):
    def __init__(self, symbols: Sequence[str], name: str) -> None:
        self.symbols: Final = symbols
        self.name: Final = name

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if doc in self.symbols:
            return doc
        raise ValidationException(f"Expected one of {self.symbols}")

    def __repr__(self) -> str:
        return self.name


class _SecondaryDSLLoader(_Loader):
    def __init__(self, inner: _Loader) -> None:
        self.inner: Final = inner

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        r: Final[list[dict[str, Any]]] = []
        if isinstance(doc, MutableSequence):
            for d in doc:
                if isinstance(d, str):
                    if d.endswith("?"):
                        r.append({"pattern": d[:-1], "required": False})
                    else:
                        r.append({"pattern": d})
                elif isinstance(d, dict):
                    new_dict1: dict[str, Any] = {}
                    dict_copy = copy.deepcopy(d)
                    if "pattern" in dict_copy:
                        new_dict1["pattern"] = dict_copy.pop("pattern")
                    else:
                        raise ValidationException(
                            f"Missing pattern in secondaryFiles specification entry: {d}"
                        )
                    new_dict1["required"] = (
                        dict_copy.pop("required") if "required" in dict_copy else None
                    )

                    if len(dict_copy):
                        raise ValidationException(
                            "Unallowed values in secondaryFiles specification entry: {}".format(
                                dict_copy
                            )
                        )
                    r.append(new_dict1)

                else:
                    raise ValidationException(
                        "Expected a string or sequence of (strings or mappings)."
                    )
        elif isinstance(doc, MutableMapping):
            new_dict2: Final = {}
            doc_copy: Final = copy.deepcopy(doc)
            if "pattern" in doc_copy:
                new_dict2["pattern"] = doc_copy.pop("pattern")
            else:
                raise ValidationException(
                    f"Missing pattern in secondaryFiles specification entry: {doc}"
                )
            new_dict2["required"] = doc_copy.pop("required") if "required" in doc_copy else None

            if len(doc_copy):
                raise ValidationException(
                    f"Unallowed values in secondaryFiles specification entry: {doc_copy}"
                )
            r.append(new_dict2)

        elif isinstance(doc, str):
            if doc.endswith("?"):
                r.append({"pattern": doc[:-1], "required": False})
            else:
                r.append({"pattern": doc})
        else:
            raise ValidationException("Expected str or sequence of str")
        return self.inner.load(r, baseuri, loadingOptions, docRoot, lc=lc)


class _RecordLoader(_Loader):
    def __init__(
        self,
        classtype: type[Saveable],
        container: Optional[str] = None,
        no_link_check: Optional[bool] = None,
    ) -> None:
        self.classtype: Final = classtype
        self.container: Final = container
        self.no_link_check: Final = no_link_check

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if not isinstance(doc, MutableMapping):
            raise ValidationException(
                f"Value is a {convert_typing(extract_type(type(doc)))}, "
                f"but valid type for this field is an object."
            )
        if self.container is not None or self.no_link_check is not None:
            loadingOptions = LoadingOptions(
                copyfrom=loadingOptions, container=self.container, no_link_check=self.no_link_check
            )
        return self.classtype.fromDoc(doc, baseuri, loadingOptions, docRoot=docRoot)

    def __repr__(self) -> str:
        return str(self.classtype.__name__)


class _ExpressionLoader(_Loader):
    def __init__(self, items: type[str]) -> None:
        self.items: Final = items

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if not isinstance(doc, str):
            raise ValidationException(
                f"Value is a {convert_typing(extract_type(type(doc)))}, "
                f"but valid type for this field is a str."
            )
        return doc


class _UnionLoader(_Loader):
    def __init__(self, alternates: Sequence[_Loader], name: Optional[str] = None) -> None:
        self.alternates = alternates
        self.name: Final = name

    def add_loaders(self, loaders: Sequence[_Loader]) -> None:
        self.alternates = tuple(loader for loader in chain(self.alternates, loaders))

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        errors: Final = []

        if lc is None:
            lc = []

        for t in self.alternates:
            try:
                return t.load(doc, baseuri, loadingOptions, docRoot=docRoot, lc=lc)
            except ValidationException as e:
                if isinstance(t, _ArrayLoader) and len(self.alternates) > 1:
                    continue
                if isinstance(doc, (CommentedMap, dict)):
                    if "class" in doc:
                        if str(doc.get("class")) == str(t):
                            errors.append(
                                ValidationException(
                                    f"Object `{baseuri.split('/')[-1]}` is not valid because:",
                                    SourceLine(doc, next(iter(doc)), str),
                                    [e],
                                )
                            )
                    else:
                        if "array" in str(t):
                            continue
                        else:
                            if "id" in doc:
                                id = baseuri.split("/")[-1] + "#" + str(doc.get("id"))
                                if "id" in lc:
                                    errors.append(
                                        ValidationException(
                                            f"checking object `{id}` using `{t}`",
                                            SourceLine(lc, "id", str),
                                            [e],
                                        )
                                    )
                                else:
                                    errors.append(
                                        ValidationException(
                                            f"checking object `{id}` using `{t}`",
                                            SourceLine(lc, doc.get("id"), str),
                                            [e],
                                        )
                                    )
                            else:
                                if not isinstance(
                                    t, (_PrimitiveLoader)
                                ):  # avoids 'tried <class "NoneType"> was {x}' errors
                                    errors.append(
                                        ValidationException(f"tried `{t}` but", None, [e])
                                    )
                else:
                    # avoids "tried <class "CWLType"> but x" and instead returns the values for parsing
                    errors.append(ValidationException("", None, [e]))

        if isinstance(doc, (CommentedMap, dict)) and "class" in doc:
            if str(doc.get("class")) not in str(self.alternates):
                errors.append(
                    ValidationException(
                        "Field `class` contains undefined reference to "
                        + "`"
                        + "/".join(baseuri.split("/")[0:-1])
                        + "/"
                        + str(doc.get("class"))
                        + "`",
                        SourceLine(doc, "class", str),
                        [],
                    )
                )
        raise ValidationException("", None, errors, "*")

    def __repr__(self) -> str:
        return self.name if self.name is not None else " | ".join(str(a) for a in self.alternates)


class _URILoader(_Loader):
    def __init__(
        self,
        inner: _Loader,
        scoped_id: bool,
        vocab_term: bool,
        scoped_ref: Optional[int],
        no_link_check: Optional[bool],
    ) -> None:
        self.inner: Final = inner
        self.scoped_id: Final = scoped_id
        self.vocab_term: Final = vocab_term
        self.scoped_ref: Final = scoped_ref
        self.no_link_check: Final = no_link_check

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if self.no_link_check is not None:
            loadingOptions = LoadingOptions(
                copyfrom=loadingOptions, no_link_check=self.no_link_check
            )
        if isinstance(doc, MutableSequence):
            newdoc: Final = []
            for i in doc:
                if isinstance(i, str):
                    newdoc.append(
                        expand_url(
                            i,
                            baseuri,
                            loadingOptions,
                            self.scoped_id,
                            self.vocab_term,
                            self.scoped_ref,
                        )
                    )
                else:
                    newdoc.append(i)
            doc = newdoc
        elif isinstance(doc, str):
            doc = expand_url(
                doc,
                baseuri,
                loadingOptions,
                self.scoped_id,
                self.vocab_term,
                self.scoped_ref,
            )
        if isinstance(doc, str):
            if not loadingOptions.no_link_check:
                errors: Final = []
                try:
                    if not loadingOptions.fetcher.check_exists(doc):
                        errors.append(
                            ValidationException(f"contains undefined reference to `{doc}`")
                        )
                except ValidationException:
                    pass
                if len(errors) > 0:
                    raise ValidationException("", None, errors)
        return self.inner.load(doc, baseuri, loadingOptions, lc=lc)


class _TypeDSLLoader(_Loader):
    def __init__(self, inner: _Loader, refScope: Optional[int], salad_version: str) -> None:
        self.inner: Final = inner
        self.refScope: Final = refScope
        self.salad_version: Final = salad_version

    def resolve(
        self,
        doc: str,
        baseuri: str,
        loadingOptions: LoadingOptions,
    ) -> Union[list[Union[dict[str, Any], str]], dict[str, Any], str]:
        doc_ = doc
        optional = False
        if doc_.endswith("?"):
            optional = True
            doc_ = doc_[0:-1]

        if doc_.endswith("[]"):
            salad_versions: Final = [int(v) for v in self.salad_version[1:].split(".")]
            items: Union[list[Union[dict[str, Any], str]], dict[str, Any], str] = ""
            rest: Final = doc_[0:-2]
            if salad_versions < [1, 3]:
                if rest.endswith("[]"):
                    # To show the error message with the original type
                    return doc
                else:
                    items = expand_url(rest, baseuri, loadingOptions, False, True, self.refScope)
            else:
                items = self.resolve(rest, baseuri, loadingOptions)
                if isinstance(items, str):
                    items = expand_url(items, baseuri, loadingOptions, False, True, self.refScope)
            expanded: Union[dict[str, Any], str] = {"type": "array", "items": items}
        else:
            expanded = expand_url(doc_, baseuri, loadingOptions, False, True, self.refScope)

        if optional:
            return ["null", expanded]
        else:
            return expanded

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if isinstance(doc, MutableSequence):
            r: Final[list[Any]] = []
            for d in doc:
                if isinstance(d, str):
                    resolved = self.resolve(d, baseuri, loadingOptions)
                    if isinstance(resolved, MutableSequence):
                        for i in resolved:
                            if i not in r:
                                r.append(i)
                    else:
                        if resolved not in r:
                            r.append(resolved)
                else:
                    r.append(d)
            doc = r
        elif isinstance(doc, str):
            doc = self.resolve(doc, baseuri, loadingOptions)

        return self.inner.load(doc, baseuri, loadingOptions, lc=lc)


class _IdMapLoader(_Loader):
    def __init__(self, inner: _Loader, mapSubject: str, mapPredicate: Optional[str]) -> None:
        self.inner: Final = inner
        self.mapSubject: Final = mapSubject
        self.mapPredicate: Final = mapPredicate

    def load(
        self,
        doc: Any,
        baseuri: str,
        loadingOptions: LoadingOptions,
        docRoot: Optional[str] = None,
        lc: Optional[list[Any]] = None,
    ) -> Any:
        if isinstance(doc, MutableMapping):
            r: Final[list[Any]] = []
            for k in doc.keys():
                val = doc[k]
                if isinstance(val, CommentedMap):
                    v = copy.copy(val)
                    v.lc.data = val.lc.data
                    v.lc.filename = val.lc.filename
                    v[self.mapSubject] = k
                    r.append(v)
                elif isinstance(val, MutableMapping):
                    v2 = copy.copy(val)
                    v2[self.mapSubject] = k
                    r.append(v2)
                else:
                    if self.mapPredicate:
                        v3 = {self.mapPredicate: val}
                        v3[self.mapSubject] = k
                        r.append(v3)
                    else:
                        raise ValidationException("No mapPredicate")
            doc = r
        return self.inner.load(doc, baseuri, loadingOptions, lc=lc)


def _document_load(
    loader: _Loader,
    doc: Union[str, MutableMapping[str, Any], MutableSequence[Any]],
    baseuri: str,
    loadingOptions: LoadingOptions,
    addl_metadata_fields: Optional[MutableSequence[str]] = None,
) -> tuple[Any, LoadingOptions]:
    if isinstance(doc, str):
        return _document_load_by_url(
            loader,
            loadingOptions.fetcher.urljoin(baseuri, doc),
            loadingOptions,
            addl_metadata_fields=addl_metadata_fields,
        )

    if isinstance(doc, MutableMapping):
        addl_metadata: Final = {}
        if addl_metadata_fields is not None:
            for mf in addl_metadata_fields:
                if mf in doc:
                    addl_metadata[mf] = doc[mf]

        docuri: Final = baseuri
        if "$base" in doc:
            baseuri = doc["$base"]

        loadingOptions = LoadingOptions(
            copyfrom=loadingOptions,
            namespaces=doc.get("$namespaces", None),
            schemas=doc.get("$schemas", None),
            baseuri=doc.get("$base", None),
            addl_metadata=addl_metadata,
        )

        doc2: Final = copy.copy(doc)
        if "$namespaces" in doc2:
            doc2.pop("$namespaces")
        if "$schemas" in doc2:
            doc2.pop("$schemas")
        if "$base" in doc2:
            doc2.pop("$base")

        if "$graph" in doc2:
            loadingOptions.idx[baseuri] = (
                loader.load(doc2["$graph"], baseuri, loadingOptions),
                loadingOptions,
            )
        else:
            loadingOptions.idx[baseuri] = (
                loader.load(doc2, baseuri, loadingOptions, docRoot=baseuri),
                loadingOptions,
            )

        if docuri != baseuri:
            loadingOptions.idx[docuri] = loadingOptions.idx[baseuri]

        return loadingOptions.idx[baseuri]

    if isinstance(doc, MutableSequence):
        loadingOptions.idx[baseuri] = (
            loader.load(doc, baseuri, loadingOptions),
            loadingOptions,
        )
        return loadingOptions.idx[baseuri]

    raise ValidationException(
        "Expected URI string, MutableMapping or MutableSequence, got %s" % type(doc)
    )


def _document_load_by_url(
    loader: _Loader,
    url: str,
    loadingOptions: LoadingOptions,
    addl_metadata_fields: Optional[MutableSequence[str]] = None,
) -> tuple[Any, LoadingOptions]:
    if url in loadingOptions.idx:
        return loadingOptions.idx[url]

    doc_url, frg = urldefrag(url)

    text: Final = loadingOptions.fetcher.fetch_text(doc_url)
    textIO: Final = StringIO(text)
    textIO.name = str(doc_url)
    yaml: Final = yaml_no_ts()
    result: Final = yaml.load(textIO)
    add_lc_filename(result, doc_url)

    loadingOptions = LoadingOptions(copyfrom=loadingOptions, fileuri=doc_url)

    _document_load(
        loader,
        result,
        doc_url,
        loadingOptions,
        addl_metadata_fields=addl_metadata_fields,
    )

    return loadingOptions.idx[url]


def file_uri(path: str, split_frag: bool = False) -> str:
    """Transform a file path into a URL with file scheme."""
    if path.startswith("file://"):
        return path
    if split_frag:
        pathsp: Final = path.split("#", 2)
        frag = "#" + quote(str(pathsp[1])) if len(pathsp) == 2 else ""
        urlpath = pathname2url(str(pathsp[0]))
    else:
        urlpath = pathname2url(path)
        frag = ""
    if urlpath.startswith("//"):
        return f"file:{urlpath}{frag}"
    return f"file://{urlpath}{frag}"


def prefix_url(url: str, namespaces: dict[str, str]) -> str:
    """Expand short forms into full URLs using the given namespace dictionary."""
    for k, v in namespaces.items():
        if url.startswith(v):
            return k + ":" + url[len(v) :]
    return url


def save_relative_uri(
    uri: Any,
    base_url: str,
    scoped_id: bool,
    ref_scope: Optional[int],
    relative_uris: bool,
) -> Any:
    """Convert any URI to a relative one, obeying the scoping rules."""
    if isinstance(uri, MutableSequence):
        return [save_relative_uri(u, base_url, scoped_id, ref_scope, relative_uris) for u in uri]
    elif isinstance(uri, str):
        if not relative_uris or uri == base_url:
            return uri
        urisplit: Final = urlsplit(uri)
        basesplit: Final = urlsplit(base_url)
        if urisplit.scheme == basesplit.scheme and urisplit.netloc == basesplit.netloc:
            if urisplit.path != basesplit.path:
                p = os.path.relpath(urisplit.path, os.path.dirname(basesplit.path))
                if urisplit.fragment:
                    p = p + "#" + urisplit.fragment
                return p

            basefrag = basesplit.fragment + "/"
            if ref_scope:
                sp = basefrag.split("/")
                i = 0
                while i < ref_scope:
                    sp.pop()
                    i += 1
                basefrag = "/".join(sp)

            if urisplit.fragment.startswith(basefrag):
                return urisplit.fragment[len(basefrag) :]
            return urisplit.fragment
        return uri
    else:
        return save(uri, top=False, base_url=base_url, relative_uris=relative_uris)


def shortname(inputid: str) -> str:
    """
    Compute the shortname of a fully qualified identifier.

    See https://w3id.org/cwl/v1.2/SchemaSalad.html#Short_names.
    """
    parsed_id: Final = urlparse(inputid)
    if parsed_id.fragment:
        return parsed_id.fragment.split("/")[-1]
    return parsed_id.path.split("/")[-1]
