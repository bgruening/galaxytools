import copy
import logging
import os
import pathlib
import re
import tempfile
import traceback
import urllib
import xml.sax  # nosec
from collections.abc import MutableMapping, MutableSequence
from io import StringIO
from typing import Any, Callable, Final, Optional, Union, cast

import requests
from cachecontrol.caches import SeparateBodyFileCache
from cachecontrol.wrapper import CacheControl
from mypy_extensions import mypyc_attr
from rdflib.exceptions import ParserError
from rdflib.graph import Graph
from rdflib.namespace import OWL, RDF, RDFS
from rdflib.plugin import PluginException
from rdflib.plugins.parsers.notation3 import BadSyntax
from rdflib.util import guess_format
from ruamel.yaml.comments import CommentedMap, CommentedSeq, LineCol
from ruamel.yaml.error import MarkedYAMLError

from .exceptions import SchemaSaladException, ValidationException
from .fetcher import DefaultFetcher
from .sourceline import SourceLine, add_lc_filename, relname
from .utils import (
    AttachmentsType,
    CacheType,
    ContextType,
    FetcherCallableType,
    IdxResultType,
    IdxType,
    ResolvedRefType,
    ResolveType,
    aslist,
    onWindows,
    yaml_no_ts,
)

_logger: Final = logging.getLogger("salad")


def file_uri(path: str, split_frag: bool = False) -> str:
    if path.startswith("file://"):
        return path
    if split_frag:
        pathsp: Final = path.split("#", 2)
        if len(pathsp) == 2:
            frag = "#" + urllib.parse.quote(str(pathsp[1]))
        else:
            frag = ""
        urlpath = urllib.request.pathname2url(str(pathsp[0]))
    else:
        urlpath = urllib.request.pathname2url(path)
        frag = ""
    if urlpath.startswith("//"):
        return f"file:{urlpath}{frag}"
    return f"file://{urlpath}{frag}"


def uri_file_path(url: str) -> str:
    split: Final = urllib.parse.urlsplit(url)
    if split.scheme == "file":
        return urllib.request.url2pathname(str(split.path)) + (
            "#" + urllib.parse.unquote(str(split.fragment)) if bool(split.fragment) else ""
        )
    raise ValidationException(f"Not a file URI: {url}")


def to_validation_exception(e: MarkedYAMLError) -> ValidationException:
    """Convert ruamel.yaml exception to our type."""
    fname_regex: Final = re.compile(r"^file://" + re.escape(os.getcwd()) + "/")

    exc: Final = ValidationException(e.problem)
    mark: Final = e.problem_mark
    exc.file = re.sub(fname_regex, "", mark.name)
    exc.start = (mark.line + 1, mark.column + 1)
    exc.end = None

    if e.context:
        parent: Final = ValidationException(e.context)
        context_mark: Final = e.context_mark
        if context_mark:
            parent.file = re.sub(fname_regex, "", context_mark.name)
            parent.start = (context_mark.line + 1, context_mark.column + 1)
        parent.end = None
        parent.children = [exc]
        return parent
    return exc


class NormDict(dict[str, Union[CommentedMap, CommentedSeq, str, None]]):
    """A Dict where all keys are normalized using the provided function."""

    def __init__(self, normalize: Callable[[str], str] = str) -> None:
        super().__init__()
        self.normalize: Final = normalize

    def __eq__(self, other: Any) -> bool:
        return super().__eq__(other)

    def __getitem__(self, key: Any) -> Any:
        return super().__getitem__(self.normalize(key))

    def __setitem__(self, key: Any, value: Any) -> Any:
        return super().__setitem__(self.normalize(key), value)

    def __delitem__(self, key: Any) -> Any:
        return super().__delitem__(self.normalize(key))

    def __contains__(self, key: Any) -> bool:
        return super().__contains__(self.normalize(key))

    def __del__(self) -> None:
        del self.normalize


def SubLoader(loader: "Loader") -> "Loader":
    return Loader(
        loader.ctx,
        schemagraph=loader.graph,
        foreign_properties=loader.foreign_properties,
        idx=loader.idx,
        cache=loader.cache,
        fetcher_constructor=loader.fetcher_constructor,
        skip_schemas=loader.skip_schemas,
        url_fields=loader.url_fields,
        allow_attachments=loader.allow_attachments,
        session=loader.session,
        salad_version=loader.salad_version,
    )


def _url_norm(url: str) -> str:
    return urllib.parse.urlsplit(url).geturl()


@mypyc_attr(allow_interpreted_subclasses=True)
class Loader:
    def __init__(
        self,
        ctx: ContextType,
        schemagraph: Optional[Graph] = None,
        foreign_properties: Optional[set[str]] = None,
        idx: Optional[IdxType] = None,
        cache: Optional[CacheType] = None,
        session: Optional[requests.sessions.Session] = None,
        fetcher_constructor: Optional[FetcherCallableType] = None,
        skip_schemas: Optional[bool] = None,
        url_fields: Optional[set[str]] = None,
        allow_attachments: Optional[AttachmentsType] = None,
        doc_cache: Union[str, bool] = True,
        salad_version: Optional[str] = None,
    ) -> None:
        self.idx: IdxType = NormDict(_url_norm) if idx is None else idx

        self.ctx: Final[ContextType] = {}
        self.graph = schemagraph if schemagraph is not None else Graph()
        self.foreign_properties: Final = (
            set(foreign_properties) if foreign_properties is not None else set()
        )
        self.cache: Final = cache if cache is not None else {}
        self.skip_schemas: Final = skip_schemas if skip_schemas is not None else False

        if session is None:
            if doc_cache is False:
                temp_session = requests.Session()
            elif doc_cache is True:
                root: Final = pathlib.Path(os.environ.get("HOME", tempfile.gettempdir()))
                temp_session = CacheControl(
                    requests.Session(),
                    cache=SeparateBodyFileCache(root / ".cache" / "salad"),
                )
            elif isinstance(doc_cache, str):
                temp_session = CacheControl(
                    requests.Session(), cache=SeparateBodyFileCache(doc_cache)
                )
        else:
            temp_session = session
        self.session = temp_session

        self.fetcher_constructor: Final = (
            fetcher_constructor if fetcher_constructor is not None else DefaultFetcher
        )
        self.fetcher: Final = self.fetcher_constructor(self.cache, self.session)
        self.fetch_text: Final = self.fetcher.fetch_text
        self.check_exists: Final = self.fetcher.check_exists
        self.url_fields: Final[set[str]] = set() if url_fields is None else set(url_fields)
        self.scoped_ref_fields: Final[dict[str, int]] = {}
        self.vocab_fields: Final[set[str]] = set()
        self.identifiers: Final[list[str]] = []
        self.identity_links: Final[set[str]] = set()
        self.standalone: Optional[set[str]] = None
        self.nolinkcheck: Final[set[str]] = set()
        self.vocab: Final[dict[str, str]] = {}
        self.rvocab: Final[dict[str, str]] = {}
        self.idmap: Final[dict[str, str]] = {}
        self.mapPredicate: Final[dict[str, str]] = {}
        self.type_dsl_fields: Final[set[str]] = set()
        self.subscopes: Final[dict[str, str]] = {}
        self.secondaryFile_dsl_fields: Final[set[str]] = set()
        self.allow_attachments: Final = allow_attachments

        if salad_version:
            temp_salad_version = salad_version
        else:
            temp_salad_version = "v1.1"
        self.salad_version: Final = temp_salad_version

        self.add_context(ctx)

    def expand_url(
        self,
        url: str,
        base_url: str,
        scoped_id: bool = False,
        vocab_term: bool = False,
        scoped_ref: Optional[int] = None,
    ) -> str:
        if url in ("@id", "@type"):
            return url

        if vocab_term and url in self.vocab:
            return url

        if url.startswith("_:"):
            return url

        if bool(self.vocab) and ":" in url:
            prefix: Final = url.split(":")[0]
            if not prefix:
                pass
            elif prefix in self.vocab:
                url = self.vocab[prefix] + url[len(prefix) + 1 :]
            elif (
                prefix not in self.fetcher.supported_schemes()
                and "/" not in prefix
                and "#" not in prefix
            ):
                _logger.warning(
                    "URI prefix '%s' of '%s' not recognized, are you missing a "
                    "$namespaces section?",
                    prefix,
                    url,
                )

        split: Final = urllib.parse.urlsplit(url)

        if (
            (bool(split.scheme) and split.scheme in self.fetcher.supported_schemes())
            or url.startswith("$(")
            or url.startswith("${")
        ):
            pass
        elif scoped_id and not bool(split.fragment):
            splitbase: Final = urllib.parse.urlsplit(base_url)
            frg: Final = (
                splitbase.fragment + "/" + split.path if bool(splitbase.fragment) else split.path
            )
            pt: Final = splitbase.path if splitbase.path != "" else "/"
            url = urllib.parse.urlunsplit(
                (splitbase.scheme, splitbase.netloc, pt, splitbase.query, frg)
            )
        elif scoped_ref is not None and not split.fragment:
            pass
        else:
            url = self.fetcher.urljoin(base_url, url)

        if vocab_term and url in self.rvocab:
            return self.rvocab[url]
        return url

    def _add_properties(self, s: str) -> None:
        for _, _, rng in self.graph.triples((s, RDFS.range, None)):
            literal = (
                str(rng).startswith("http://www.w3.org/2001/XMLSchema#")
                and str(rng) != "http://www.w3.org/2001/XMLSchema#anyURI"
            ) or str(rng) == "http://www.w3.org/2000/01/rdf-schema#Literal"
            if not literal:
                self.url_fields.add(str(s))
        self.foreign_properties.add(str(s))

    def add_namespaces(self, ns: dict[str, str]) -> None:
        """Add the given namespace to our vocab list."""
        self.vocab.update(ns)

    def add_schemas(self, ns: Union[list[str], str], base_url: str) -> None:
        """Fetch external schemas and add them to the graph."""
        if self.skip_schemas:
            return
        for sch in aslist(ns):
            fetchurl = self.fetcher.urljoin(base_url, sch)
            if fetchurl not in self.cache or self.cache[fetchurl] is True:
                _logger.debug("Getting external schema %s", fetchurl)
                try:
                    content = self.fetch_text(fetchurl)
                except Exception as e:
                    tb = traceback.format_exception(type(e), e, e.__traceback__)
                    _logger.warning("Could not load extension schema %s: %s", fetchurl, str(e))
                    _logger.debug(tb)
                    continue
                newGraph = Graph()
                err_msg = "unknown error"
                for fmt in [
                    guess_format(sch),
                    "xml",
                    "turtle",
                    None,
                    guess_format(sch),
                ]:
                    try:
                        newGraph.parse(data=content, format=fmt, publicID=str(fetchurl))
                        self.cache[fetchurl] = newGraph
                        self.graph += newGraph
                        break
                    except (
                        xml.sax.SAXParseException,
                        TypeError,
                        BadSyntax,
                        ParserError,
                        PluginException,
                    ) as e:
                        err_msg = str(e)
                else:
                    _logger.warning("Could not load extension schema %s: %s", fetchurl, err_msg)

        for s, _, _ in self.graph.triples((None, RDF.type, RDF.Property)):
            self._add_properties(s)
        for s, _, o in self.graph.triples((None, RDFS.subPropertyOf, None)):
            self._add_properties(s)
            self._add_properties(o)
        for s, _, _ in self.graph.triples((None, RDFS.range, None)):
            self._add_properties(s)
        for s, _, _ in self.graph.triples((None, RDF.type, OWL.ObjectProperty)):
            self._add_properties(s)

        for s, _, _ in self.graph.triples((None, None, None)):
            self.idx[str(s)] = None

    def add_context(self, newcontext: ContextType) -> None:
        if bool(self.vocab):
            raise ValidationException("Refreshing context that already has stuff in it")

        self.url_fields.add("$schemas")
        self.scoped_ref_fields.clear()
        self.vocab_fields.clear()
        self.identifiers.clear()
        self.identity_links.clear()
        self.standalone = set()
        self.nolinkcheck.clear()
        self.idmap.clear()
        self.mapPredicate.clear()
        self.vocab.clear()
        self.rvocab.clear()
        self.type_dsl_fields.clear()
        self.secondaryFile_dsl_fields.clear()
        self.subscopes.clear()

        self.ctx.update(_copy_dict_without_key(newcontext, "@context"))

        _logger.debug("ctx is %s", self.ctx)

        for key, value in self.ctx.items():
            if value == "@id":
                self.identifiers.append(key)
                self.identity_links.add(key)
            elif isinstance(value, MutableMapping):
                if value.get("@type") == "@id":
                    self.url_fields.add(key)
                    if "refScope" in value:
                        self.scoped_ref_fields[key] = value["refScope"]
                    if value.get("identity", False):
                        self.identity_links.add(key)

                if value.get("@type") == "@vocab":
                    self.url_fields.add(key)
                    self.vocab_fields.add(key)
                    if "refScope" in value:
                        self.scoped_ref_fields[key] = value["refScope"]
                    if value.get("typeDSL"):
                        self.type_dsl_fields.add(key)

                if value.get("secondaryFilesDSL"):
                    self.secondaryFile_dsl_fields.add(key)

                if value.get("noLinkCheck"):
                    self.nolinkcheck.add(key)

                if value.get("mapSubject"):
                    self.idmap[key] = value["mapSubject"]

                if value.get("mapPredicate"):
                    self.mapPredicate[key] = value["mapPredicate"]

                if value.get("@id"):
                    self.vocab[key] = value["@id"]

                if value.get("subscope"):
                    self.subscopes[key] = value["subscope"]

            elif isinstance(value, str):
                self.vocab[key] = value

        for k, v in self.vocab.items():
            self.rvocab[self.expand_url(v, "", scoped_id=False)] = k

        self.identifiers.sort()

        _logger.debug("identifiers is %s", self.identifiers)
        _logger.debug("identity_links is %s", self.identity_links)
        _logger.debug("url_fields is %s", self.url_fields)
        _logger.debug("vocab_fields is %s", self.vocab_fields)
        _logger.debug("vocab is %s", self.vocab)

    def resolve_ref(
        self,
        ref: ResolveType,
        base_url: Optional[str] = None,
        checklinks: bool = True,
        strict_foreign_properties: bool = False,
        content_types: Optional[list[str]] = None,  # Expected content-types
    ) -> ResolvedRefType:
        lref = ref
        obj: Optional[CommentedMap] = None
        resolved_obj: ResolveType = None
        imp = False
        inc = False
        mixin: Optional[MutableMapping[str, str]] = None

        if not base_url:
            base_url = file_uri(os.getcwd()) + "/"

        # If `ref` is a dict, look for special directives.
        if isinstance(lref, CommentedMap):
            obj = lref
            if "$import" in obj:
                if len(obj) == 1:
                    lref = obj["$import"]
                    imp = True
                    obj = None
                else:
                    raise ValidationException(
                        f"'$import' must be the only field in {obj}",
                        SourceLine(obj, "$import"),
                    )
            elif "$include" in obj:
                if len(obj) == 1:
                    lref = obj["$include"]
                    inc = True
                    obj = None
                else:
                    raise ValidationException(
                        f"'$include' must be the only field in {obj}",
                        SourceLine(obj, "$include"),
                    )
            elif "$mixin" in obj:
                lref = obj["$mixin"]
                mixin = obj
                obj = None
            else:
                lref = None
                for identifier in self.identifiers:
                    if identifier in obj:
                        lref = obj[identifier]
                        break
                if not lref:
                    raise ValidationException(
                        f"Object {obj!r} does not have identifier field in {self.identifiers}",
                        SourceLine(obj),
                    )

        if not isinstance(lref, str):
            raise ValidationException(
                f"Expected CommentedMap or string, got {type(lref)}: {lref!r}"
            )

        if isinstance(lref, str) and os.sep == "\\":
            # Convert Windows path separator in ref
            lref = lref.replace("\\", "/")

        url = self.expand_url(lref, base_url, scoped_id=(obj is not None))
        # Has this reference been loaded already?
        if url in self.idx and (not mixin):
            resolved_obj = self.idx[url]
            if isinstance(resolved_obj, MutableMapping):
                metadata: Union[CommentedMap, CommentedSeq, str, None] = self.idx.get(
                    urllib.parse.urldefrag(url)[0], CommentedMap()
                )
                if isinstance(metadata, MutableMapping):
                    if "$graph" in resolved_obj:
                        metadata = _copy_dict_without_key(resolved_obj, "$graph")
                        return resolved_obj["$graph"], metadata
                    else:
                        return resolved_obj, metadata
                else:
                    raise ValidationException(
                        f"Expected CommentedMap, got {type(metadata)}: {metadata!r}"
                    )
            elif isinstance(resolved_obj, MutableSequence):
                metadata = self.idx.get(urllib.parse.urldefrag(url)[0], CommentedMap())
                if isinstance(metadata, MutableMapping):
                    return resolved_obj, metadata
                else:
                    return resolved_obj, CommentedMap()
            elif isinstance(resolved_obj, str):
                return resolved_obj, CommentedMap()
            else:
                raise ValidationException(
                    f"Expected MutableMapping or MutableSequence, got "
                    f"{type(resolved_obj)}: {resolved_obj!r}"
                )

        # "$include" directive means load raw text
        if inc:
            # Make a note in the index that this was an included string
            self.idx["include:" + url] = url
            return self.fetch_text(url), CommentedMap()

        doc = None
        if isinstance(obj, MutableMapping):
            for identifier in self.identifiers:
                obj[identifier] = url
            doc_url = url
        else:
            # Load structured document
            doc_url, frg = urllib.parse.urldefrag(url)
            if doc_url in self.idx and (not mixin):
                # If the base document is in the index, it was already loaded,
                # so if we didn't find the reference earlier then it must not
                # exist.
                raise SourceLine(self.idx, doc_url, ValidationException).makeError(
                    f"Reference '#{frg}' not found in file {doc_url!r}.",
                )
            doc = self.fetch(doc_url, inject_ids=(not mixin), content_types=content_types)

        if imp:
            # Make a note in the index that this was an imported fragment
            self.idx["import:" + url] = url

        # Recursively expand urls and resolve directives
        if bool(mixin):
            doc = copy.deepcopy(doc)
            if isinstance(doc, CommentedMap) and mixin is not None:
                doc.update(mixin)
                del doc["$mixin"]
            resolved_obj, metadata = self.resolve_all(
                doc,
                base_url,
                file_base=doc_url,
                checklinks=checklinks,
                strict_foreign_properties=strict_foreign_properties,
            )
        else:
            resolved_obj, metadata = self.resolve_all(
                doc or obj,
                doc_url,
                checklinks=checklinks,
                strict_foreign_properties=strict_foreign_properties,
            )

        # Requested reference should be in the index now, otherwise it's a bad
        # reference
        if not bool(mixin):
            if url in self.idx:
                resolved_obj = self.idx[url]
            else:
                raise ValidationException(
                    f"Reference {url!r} is not in the index. Index contains: {' '.join(self.idx)}"
                )

        if isinstance(resolved_obj, CommentedMap):
            if "$graph" in resolved_obj:
                metadata = _copy_dict_without_key(resolved_obj, "$graph")
                return resolved_obj["$graph"], metadata
            else:
                return resolved_obj, metadata
        else:
            return resolved_obj, metadata

    def _resolve_idmap(
        self,
        document: CommentedMap,
        loader: "Loader",
    ) -> None:
        # Convert fields with mapSubject into lists
        # use mapPredicate if the mapped value isn't a dict.
        for idmapField in loader.idmap:
            if idmapField in document:
                idmapFieldValue = document[idmapField]
                if (
                    isinstance(idmapFieldValue, MutableMapping)
                    and "$import" not in idmapFieldValue
                    and "$include" not in idmapFieldValue
                ):
                    ls = CommentedSeq()
                    for k in sorted(idmapFieldValue.keys()):
                        val = idmapFieldValue[k]
                        v: Optional[CommentedMap] = None
                        if not isinstance(val, CommentedMap):
                            if idmapField in loader.mapPredicate:
                                v = CommentedMap(((loader.mapPredicate[idmapField], val),))
                                v.lc.add_kv_line_col(
                                    loader.mapPredicate[idmapField],
                                    document[idmapField].lc.data[k],
                                )
                                v.lc.filename = document.lc.filename
                            else:
                                raise ValidationException(
                                    f"mapSubject {k!r} value {v!r} is not a dict "
                                    "and does not have a mapPredicate.",
                                    SourceLine(document, idmapField),
                                )
                        else:
                            v = val

                        v[loader.idmap[idmapField]] = k
                        v.lc.add_kv_line_col(
                            loader.idmap[idmapField], document[idmapField].lc.data[k]
                        )
                        v.lc.filename = document.lc.filename

                        ls.lc.add_kv_line_col(len(ls), document[idmapField].lc.data[k])

                        ls.lc.filename = document.lc.filename
                        ls.append(v)

                    document[idmapField] = ls

    def _type_dsl(
        self,
        t: Union[str, CommentedMap, CommentedSeq],
        lc: LineCol,
        filename: str,
    ) -> Union[str, CommentedMap, CommentedSeq]:
        if not isinstance(t, str):
            return t

        t_ = t
        optional = False
        if t_.endswith("?"):
            optional = True
            t_ = t_[0:-1]

        if t_.endswith("[]"):
            salad_versions = [int(v) for v in self.salad_version[1:].split(".")]
            rest: Final = t_[0:-2]
            if salad_versions < [1, 3]:
                if rest.endswith("[]"):
                    # To show the error message with the original type
                    return t
                else:
                    cmap = CommentedMap((("type", "array"), ("items", rest)))
            else:
                items = self._type_dsl(rest, lc, filename)
                cmap = CommentedMap((("type", "array"), ("items", items)))
            cmap.lc.add_kv_line_col("type", lc)
            cmap.lc.add_kv_line_col("items", lc)
            cmap.lc.filename = filename
            expanded: Union[str, CommentedMap, CommentedSeq] = cmap
        else:
            expanded = t_

        if optional:
            cs = CommentedSeq(["null", expanded])
            cs.lc.add_kv_line_col(0, lc)
            cs.lc.add_kv_line_col(1, lc)
            cs.lc.filename = filename
            ret: Union[str, CommentedMap, CommentedSeq] = cs
        else:
            ret = expanded

        return ret

    def _secondaryFile_dsl(
        self,
        t: Union[str, CommentedMap, CommentedSeq],
        lc: LineCol,
        filename: str,
    ) -> Union[str, CommentedMap, CommentedSeq]:
        if not isinstance(t, str):
            return t
        pat: Final = t[0:-1] if t.endswith("?") else t
        req: Final[Optional[bool]] = False if t.endswith("?") else None

        second: Final = CommentedMap((("pattern", pat), ("required", req)))
        second.lc.add_kv_line_col("pattern", lc)
        second.lc.add_kv_line_col("required", lc)
        second.lc.filename = filename
        return second

    def _apply_dsl(
        self,
        datum: Union[str, CommentedMap, CommentedSeq],
        d: str,
        loader: "Loader",
        lc: LineCol,
        filename: str,
    ) -> Union[str, CommentedMap, CommentedSeq]:
        if d in loader.type_dsl_fields:
            return self._type_dsl(datum, lc, filename)
        if d in loader.secondaryFile_dsl_fields:
            return self._secondaryFile_dsl(datum, lc, filename)
        return datum

    def _resolve_dsl(
        self,
        document: CommentedMap,
        loader: "Loader",
    ) -> None:
        fields: Final = list(loader.type_dsl_fields)
        fields.extend(loader.secondaryFile_dsl_fields)

        for d in fields:
            if d in document:
                datum2 = datum = document[d]
                if isinstance(datum, str):
                    datum2 = self._apply_dsl(
                        datum,
                        d,
                        loader,
                        document.lc.data[d] if document.lc.data else document.lc,
                        getattr(document.lc, "filename", ""),
                    )
                elif isinstance(datum, CommentedSeq):
                    datum2 = CommentedSeq()
                    for n, t in enumerate(datum):
                        if datum.lc and datum.lc.data:
                            datum2.lc.add_kv_line_col(len(datum2), datum.lc.data[n])
                            datum2.append(
                                self._apply_dsl(
                                    t, d, loader, datum.lc.data[n], document.lc.filename
                                )
                            )
                        else:
                            datum2.append(self._apply_dsl(t, d, loader, LineCol(), ""))
                if isinstance(datum2, CommentedSeq):
                    datum3 = CommentedSeq()
                    seen: list[str] = []
                    for i, item in enumerate(datum2):
                        if isinstance(item, CommentedSeq):
                            for j, v in enumerate(item):
                                if v not in seen:
                                    datum3.lc.add_kv_line_col(len(datum3), item.lc.data[j])
                                    datum3.append(v)
                                    seen.append(v)
                        else:
                            if item not in seen:
                                if datum2.lc and datum2.lc.data:
                                    datum3.lc.add_kv_line_col(len(datum3), datum2.lc.data[i])
                                datum3.append(item)
                                seen.append(item)
                    document[d] = datum3
                else:
                    document[d] = datum2

    def _resolve_identifier(self, document: CommentedMap, loader: "Loader", base_url: str) -> str:
        # Expand identifier field (usually 'id') to resolve scope
        for identifier in loader.identifiers:
            if identifier in document:
                if isinstance(document[identifier], str):
                    document[identifier] = loader.expand_url(
                        document[identifier], base_url, scoped_id=True
                    )
                    if document[identifier] not in loader.idx or isinstance(
                        loader.idx[document[identifier]], str
                    ):
                        loader.idx[document[identifier]] = document
                    base_url = document[identifier]
                else:
                    raise ValidationException(
                        f"identifier field {document[identifier]!r} must be a string"
                    )
        return base_url

    def _resolve_identity(
        self,
        document: dict[str, Union[str, MutableSequence[Union[str, CommentedMap]]]],
        loader: "Loader",
        base_url: str,
    ) -> None:
        # Resolve scope for identity fields (fields where the value is the
        # identity of a standalone node, such as enum symbols)
        for identifier in loader.identity_links:
            if identifier in document and isinstance(document[identifier], MutableSequence):
                for n, v in enumerate(document[identifier]):
                    if isinstance(v, str):
                        document[identifier][n] = loader.expand_url(  # type: ignore
                            v, base_url, scoped_id=True
                        )
                        if document[identifier][n] not in loader.idx:
                            loader.idx[cast(str, document[identifier][n])] = v

    def _normalize_fields(self, document: CommentedMap, loader: "Loader") -> None:
        # Normalize fields which are prefixed or full URIn to vocabulary terms
        for d in list(document.keys()):
            if isinstance(d, str):
                d2 = loader.expand_url(d, "", scoped_id=False, vocab_term=True)
                if d != d2:
                    document[d2] = document[d]
                    document.lc.add_kv_line_col(d2, document.lc.data[d])
                    del document[d]

    def _resolve_uris(
        self,
        document: dict[str, Union[str, MutableSequence[Union[str, CommentedMap]]]],
        loader: "Loader",
        base_url: str,
    ) -> None:
        # Resolve remaining URLs based on document base
        for d in loader.url_fields:
            if d in document:
                datum = document[d]
                if isinstance(datum, str):
                    document[d] = loader.expand_url(
                        datum,
                        base_url,
                        scoped_id=False,
                        vocab_term=(d in loader.vocab_fields),
                        scoped_ref=loader.scoped_ref_fields.get(d),
                    )
                elif isinstance(datum, MutableSequence):
                    for i, url in enumerate(datum):
                        if isinstance(url, str):
                            datum[i] = loader.expand_url(
                                url,
                                base_url,
                                scoped_id=False,
                                vocab_term=(d in loader.vocab_fields),
                                scoped_ref=loader.scoped_ref_fields.get(d),
                            )

    def resolve_all(
        self,
        document: ResolveType,
        base_url: str,
        file_base: Optional[str] = None,
        checklinks: bool = True,
        strict_foreign_properties: bool = False,
    ) -> ResolvedRefType:
        loader = self
        metadata = CommentedMap()
        if file_base is None:
            file_base = base_url

        if isinstance(document, CommentedMap):
            # Handle $import and $include
            if "$import" in document or "$include" in document:
                return self.resolve_ref(
                    document,
                    base_url=file_base,
                    checklinks=checklinks,
                    strict_foreign_properties=strict_foreign_properties,
                )
            if "$mixin" in document:
                return self.resolve_ref(
                    document,
                    base_url=base_url,
                    checklinks=checklinks,
                    strict_foreign_properties=strict_foreign_properties,
                )
        elif isinstance(document, CommentedSeq):
            pass
        elif isinstance(document, (list, dict)):
            raise ValidationException(
                f"Expected CommentedMap or CommentedSeq, got {type(document)}: {document!r}"
            )
        else:
            return (document, metadata)

        newctx: Optional["Loader"] = None
        if isinstance(document, CommentedMap):
            # Handle $base, $profile, $namespaces, $schemas and $graph
            if "$base" in document:
                base_url = document["$base"]

            if "$profile" in document:
                if newctx is None:
                    newctx = SubLoader(self)
                newctx.add_namespaces(document.get("$namespaces", CommentedMap()))
                newctx.add_schemas(document.get("$schemas", []), document["$profile"])

            if "$namespaces" in document:
                if newctx is None:
                    newctx = SubLoader(self)
                namespaces = document["$namespaces"]
                if isinstance(namespaces, dict):
                    newctx.add_namespaces(document["$namespaces"])
                else:
                    raise ValidationException(
                        "$namespaces must be a dictionary",
                        SourceLine(document, "$namespaces"),
                    )

            if "$schemas" in document:
                if newctx is None:
                    newctx = SubLoader(self)
                schemas = document["$schemas"]
                if isinstance(schemas, (list, str)):
                    newctx.add_schemas(schemas, file_base)
                else:
                    raise ValidationException(
                        "$schemas must be a string or a list of string",
                        SourceLine(document, "$schemas"),
                    )

            if newctx is not None:
                loader = newctx

            for identifier in loader.identity_links:
                if identifier in document:
                    if isinstance(document[identifier], str):
                        document[identifier] = loader.expand_url(
                            document[identifier], base_url, scoped_id=True
                        )
                        loader.idx[document[identifier]] = document

            metadata = document
            if "$graph" in document:
                document = document["$graph"]

        if isinstance(document, CommentedMap):
            self._normalize_fields(document, loader)
            self._resolve_idmap(document, loader)
            self._resolve_dsl(document, loader)
            base_url = self._resolve_identifier(document, loader, base_url)
            self._resolve_identity(document, loader, base_url)
            self._resolve_uris(document, loader, base_url)

            try:
                for key, val in document.items():
                    subscope: str = ""
                    if key in loader.subscopes:
                        subscope = "/" + loader.subscopes[key]
                    document[key], _ = loader.resolve_all(
                        val, base_url + subscope, file_base=file_base, checklinks=False
                    )
            except ValidationException as v:
                _logger.warning("loader is %s", id(loader), exc_info=True)
                raise ValidationException(
                    f"({id(loader)}) ({file_base}) Validation error in field {key}:",
                    None,
                    [v],
                ) from v

        elif isinstance(document, CommentedSeq):
            i = 0
            try:
                while i < len(document):
                    val = document[i]
                    if isinstance(val, CommentedMap) and ("$import" in val or "$mixin" in val):
                        l, import_metadata = loader.resolve_ref(
                            val, base_url=file_base, checklinks=False
                        )
                        metadata.setdefault("$import_metadata", {})
                        for identifier in loader.identifiers:
                            if identifier in import_metadata:
                                metadata["$import_metadata"][
                                    import_metadata[identifier]
                                ] = import_metadata
                        if isinstance(l, CommentedSeq):
                            lc = document.lc.data[i]
                            del document[i]
                            llen = len(l)
                            for j in range(len(document) + llen, i + llen, -1):
                                document.lc.data[j - 1] = document.lc.data[j - llen]
                            for item in l:
                                document.insert(i, item)
                                document.lc.data[i] = lc
                                i += 1
                        else:
                            document[i] = l
                            i += 1
                    else:
                        document[i], _ = loader.resolve_all(
                            val, base_url, file_base=file_base, checklinks=False
                        )
                        i += 1
            except ValidationException as v:
                _logger.warning("failed", exc_info=True)
                raise ValidationException(
                    f"({id(loader)}) ({file_base}) Validation error in position {i}:",
                    None,
                    [v],
                ) from v

        if checklinks:
            all_doc_ids: dict[str, str] = {}
            loader.validate_links(
                document,
                "",
                all_doc_ids,
                strict_foreign_properties=strict_foreign_properties,
            )

        return document, metadata

    def fetch(
        self,
        url: str,
        inject_ids: bool = True,
        content_types: Optional[list[str]] = None,
    ) -> IdxResultType:
        if url in self.idx:
            return self.idx[url]
        try:
            text: Final = self.fetch_text(url, content_types=content_types)
            textIO: Final = StringIO(text)
            textIO.name = str(url)
            yaml: Final = yaml_no_ts()
            attachments: Final = yaml.load_all(textIO)
            result: Final = cast(Union[CommentedSeq, CommentedMap], next(attachments))

            if self.allow_attachments is not None and self.allow_attachments(result):
                i = 1
                for a in attachments:
                    self.idx[f"{url}#attachment-{i}"] = a
                    i += 1
            add_lc_filename(result, url)
        except MarkedYAMLError as e:
            raise to_validation_exception(e) from e
        if isinstance(result, CommentedMap) and inject_ids and bool(self.identifiers):
            missing_identifier = True
            for identifier in self.identifiers:
                if identifier in result:
                    missing_identifier = False
                    self.idx[self.expand_url(result[identifier], url, scoped_id=True)] = result
            if missing_identifier:
                result[self.identifiers[0]] = url
        self.idx[url] = result
        return result

    def validate_scoped(self, field: str, link: str, docid: str) -> str:
        split: Final = urllib.parse.urlsplit(docid)
        sp = split.fragment.split("/")
        n = self.scoped_ref_fields[field]
        while n > 0 and len(sp) > 0:
            sp.pop()
            n -= 1
        tried: Final = []
        while True:
            sp.append(link)
            url = urllib.parse.urlunsplit(
                (split.scheme, split.netloc, split.path, split.query, "/".join(sp))
            )
            tried.append(url)
            if url in self.idx:
                return url
            sp.pop()
            if len(sp) == 0:
                break
            sp.pop()
        if onWindows() and link.startswith("file:"):
            link = link.lower()
        raise ValidationException(
            f"Field {field!r} references unknown identifier {link!r}, tried {' '.join(tried)}"
        )

    def validate_link(
        self,
        field: str,
        link: Union[str, CommentedSeq, CommentedMap],
        # link also can be None, but that results in
        # mypyc "error: Local variable "link" has inferred type None; add an annotation"
        docid: str,
        all_doc_ids: dict[str, str],
    ) -> Union[str, CommentedSeq, CommentedMap]:
        if field in self.nolinkcheck:
            return link
        if isinstance(link, str):
            if field in self.vocab_fields:
                if link not in self.vocab and link not in self.idx and link not in self.rvocab:
                    if field in self.scoped_ref_fields:
                        return self.validate_scoped(field, link, docid)
                    elif not self.check_exists(link):
                        raise ValidationException(
                            f"Field {field!r} contains undefined reference to {link!r}"
                        )
            elif link not in self.idx and link not in self.rvocab:
                if field in self.scoped_ref_fields:
                    return self.validate_scoped(field, link, docid)
                elif not self.check_exists(link):
                    raise ValidationException(
                        f"Field {field!r} contains undefined reference to {link!r}"
                    )
        elif isinstance(link, CommentedSeq):
            errors: Final = []
            for n, i in enumerate(link):
                try:
                    link[n] = self.validate_link(field, i, docid, all_doc_ids)
                except ValidationException as v:
                    errors.append(v)
            if bool(errors):
                raise ValidationException("", None, errors)
        elif isinstance(link, CommentedMap):
            self.validate_links(link, docid, all_doc_ids)
        elif link is None:  # type: ignore[unreachable]
            return None
        else:
            raise ValidationException(
                f"{field!r} field is {type(link).__name__}, expected string, list, or a dict."
            )
        return link

    def getid(self, d: Any) -> Optional[str]:
        """Use our identifiers to extract the first match from the document."""
        if isinstance(d, MutableMapping):
            for i in self.identifiers:
                if i in d:
                    idd = d[i]
                    if isinstance(idd, str):
                        return idd
        return None

    def validate_links(
        self,
        document: ResolveType,
        base_url: str,
        all_doc_ids: dict[str, str],
        strict_foreign_properties: bool = False,
    ) -> None:
        docid: Final = self.getid(document) or base_url

        errors: Final[list[SchemaSaladException]] = []
        iterator: Any = None
        if isinstance(document, MutableSequence):
            iterator = enumerate(document)
        elif isinstance(document, MutableMapping):
            for d in self.url_fields:
                try:
                    if d in document and d not in self.identity_links:
                        document[d] = self.validate_link(d, document[d], docid, all_doc_ids)
                except SchemaSaladException as v:
                    v = v.with_sourceline(SourceLine(document, d, str))
                    if d == "$schemas" or (
                        d in self.foreign_properties and not strict_foreign_properties
                    ):
                        _logger.warning(v.as_warning())
                    else:
                        errors.append(v)
            # TODO: Validator should local scope only in which
            # duplicated keys are prohibited.
            # See also https://github.com/common-workflow-language/common-workflow-language/issues/734  # noqa: B950
            # In the future, it should raise
            # ValidationException instead of _logger.warn
            try:
                for identifier in self.identifiers:  # validate that each id is defined uniquely
                    if identifier in document:
                        sl = SourceLine(document, identifier, str)
                        if (
                            document[identifier] in all_doc_ids
                            and sl.makeLead() != all_doc_ids[document[identifier]]
                        ):
                            _logger.warning(
                                "%s object %s %r previously defined",
                                all_doc_ids[document[identifier]],
                                identifier,
                                relname(document[identifier]),
                            )
                        else:
                            all_doc_ids[document[identifier]] = sl.makeLead()
                            break
            except ValidationException as v:
                errors.append(v.with_sourceline(sl))

            iterator = list(document.items())
        else:
            return

        for key, val in iterator:
            sl = SourceLine(document, key, str)
            try:
                self.validate_links(
                    val,
                    docid,
                    all_doc_ids,
                    strict_foreign_properties=strict_foreign_properties,
                )
            except ValidationException as v:
                if key in self.nolinkcheck or (isinstance(key, str) and ":" in key):
                    _logger.warning(v.as_warning())
                else:
                    docid2 = self.getid(val)
                    if docid2 is not None:
                        errors.append(
                            ValidationException(f"checking object {relname(docid2)!r}", sl, [v])
                        )
                    else:
                        if isinstance(key, str):
                            errors.append(ValidationException(f"checking field {key!r}", sl, [v]))
                        else:
                            errors.append(ValidationException("checking item", sl, [v]))
        if bool(errors):
            if len(errors) > 1:
                raise ValidationException("", None, errors)
            raise errors[0]
        return


def _copy_dict_without_key(
    from_dict: Union[CommentedMap, ContextType], filtered_key: str
) -> CommentedMap:
    new_dict = CommentedMap(from_dict.items())
    if filtered_key in new_dict:
        del new_dict[filtered_key]
    if isinstance(from_dict, CommentedMap):
        new_dict.lc.data = copy.copy(from_dict.lc.data)
        new_dict.lc.filename = from_dict.lc.filename
    return new_dict
