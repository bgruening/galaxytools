import logging
import unicodedata
from collections.abc import Iterable, MutableMapping, MutableSequence
from typing import Any, Final, Optional, Union, cast
from urllib.parse import urldefrag, urlsplit

import rdflib
import rdflib.namespace
from rdflib import Graph, URIRef
from rdflib.namespace import RDF, RDFS
from ruamel.yaml.comments import CommentedMap, CommentedSeq

from .exceptions import SchemaException
from .utils import ContextType, aslist, json_dumps

_logger = logging.getLogger("salad")


def pred(
    datatype: MutableMapping[str, Union[dict[str, str], str]],
    field: Optional[dict[str, Any]],
    name: str,
    context: ContextType,
    defaultBase: str,
    namespaces: dict[str, rdflib.namespace.Namespace],
) -> Union[dict[str, Optional[str]], str]:
    split: Final = urlsplit(name)

    vee: Optional[str] = None

    if split.scheme != "":
        vee = name
        (ns, ln) = rdflib.namespace.split_uri(str(vee))
        name = ln
        if ns[0:-1] in namespaces:
            vee = str(namespaces[ns[0:-1]][ln])
        _logger.debug("name, v %s %s", name, vee)

    v: Optional[Union[dict[str, Optional[str]], str]] = None

    if field is not None and "jsonldPredicate" in field:
        if isinstance(field["jsonldPredicate"], MutableMapping):
            v = {}
            for k, val in field["jsonldPredicate"].items():
                v[("@" + k[1:] if k.startswith("_") else k)] = val
            if "@id" not in v:
                v["@id"] = vee
        else:
            v = field["jsonldPredicate"]
    elif "jsonldPredicate" in datatype:
        if isinstance(datatype["jsonldPredicate"], Iterable):
            for d in datatype["jsonldPredicate"]:
                if isinstance(d, MutableMapping):
                    if d["symbol"] == name:
                        v = d["predicate"]
                else:
                    raise SchemaException(
                        "entries in the jsonldPredicate List must be " "Dictionaries"
                    )
        else:
            raise SchemaException("jsonldPredicate must be a List of Dictionaries.")

    ret = v or vee

    if not ret:
        ret = defaultBase + name

    if name in context:
        if context[name] != ret:
            raise SchemaException(f"Predicate collision on {name}, {context[name]!r} != {ret!r}")
    else:
        _logger.debug("Adding to context '%s' %s (%s)", name, ret, type(ret))
        context[name] = ret

    return ret


def process_type(
    t: MutableMapping[str, Any],
    g: Graph,
    context: ContextType,
    defaultBase: str,
    namespaces: dict[str, rdflib.namespace.Namespace],
    defaultPrefix: str,
) -> None:
    if t["type"] not in ("record", "enum"):
        return

    if "name" in t:
        recordname = t["name"]

        _logger.debug("Processing %s %s\n", t.get("type"), t)

        classnode: Final = URIRef(recordname)
        g.add((classnode, RDF.type, RDFS.Class))

        split: Final = urlsplit(recordname)
        predicate = recordname
        if t.get("inVocab", True):
            if split.scheme:
                (ns, ln) = rdflib.namespace.split_uri(str(recordname))
                predicate = recordname
                recordname = ln
            else:
                predicate = f"{defaultPrefix}:{recordname}"

        if context.get(recordname, predicate) != predicate:
            raise SchemaException(
                f"Predicate collision on {recordname!r}, "
                f"{context[recordname]!r} != {predicate!r}"
            )

        if not recordname:
            raise SchemaException(f"Unable to find/derive recordname for {t}")

        _logger.debug("Adding to context '%s' %s (%s)", recordname, predicate, type(predicate))
        context[recordname] = predicate

    if t["type"] == "record":
        for i in t.get("fields", []):
            fieldname = i["name"]

            _logger.debug("Processing field %s", i)

            v: Union[dict[Any, Any], str, None] = pred(
                t, i, fieldname, context, defaultPrefix, namespaces
            )

            if isinstance(v, str):
                v = v if v[0] != "@" else None
            elif v is not None:
                v = v["_@id"] if v.get("_@id", "@")[0] != "@" else None

            if bool(v):
                try:
                    (ns, ln) = rdflib.namespace.split_uri(str(v))
                except ValueError:
                    # rdflib 5.0.0 compatibility
                    uri = str(v)
                    colon_index = str(v).rfind(":")

                    if colon_index < 0:
                        raise
                    split_start = rdflib.namespace.SPLIT_START_CATEGORIES
                    for j in range(-1 - colon_index, len(uri)):
                        if unicodedata.category(uri[j]) in split_start or uri[j] == "_":
                            # _ prevents early split, roundtrip not generate
                            ns = uri[:j]
                            if not ns:
                                break
                            ln = uri[j:]
                            break
                    if not ns or not ln:
                        raise

                if ns[0:-1] in namespaces:
                    propnode = namespaces[ns[0:-1]][ln]
                else:
                    propnode = URIRef(v)

                g.add((propnode, RDF.type, RDF.Property))
                g.add((propnode, RDFS.domain, classnode))

                # TODO generate range from datatype.

            if isinstance(i["type"], MutableMapping):
                process_type(i["type"], g, context, defaultBase, namespaces, defaultPrefix)

        if "extends" in t:
            for e in aslist(t["extends"]):
                g.add((classnode, RDFS.subClassOf, URIRef(e)))
    elif t["type"] == "enum":
        _logger.debug("Processing enum %s", t.get("name"))

        for i in t["symbols"]:
            pred(t, None, i, context, defaultBase, namespaces)


def salad_to_jsonld_context(
    j: Iterable[MutableMapping[str, Any]], schema_ctx: MutableMapping[str, Any]
) -> tuple[ContextType, Graph]:
    context: Final[ContextType] = {}
    namespaces: Final = {}
    g: Final = Graph()
    defaultPrefix: Final = ""

    for k, v in schema_ctx.items():
        context[k] = v
        namespaces[k] = rdflib.namespace.Namespace(v)

    if "@base" in context:
        defaultBase = cast(str, context["@base"])
        del context["@base"]
    else:
        defaultBase = ""

    for k, v in namespaces.items():
        g.bind(str(k), v)

    for t in j:
        process_type(t, g, context, defaultBase, namespaces, defaultPrefix)

    return (context, g)


def fix_jsonld_ids(obj: Union[CommentedMap, float, str, CommentedSeq], ids: list[str]) -> None:
    """Add missing identity entries."""
    if isinstance(obj, MutableMapping):
        for i in ids:
            if i in obj:
                obj["@id"] = obj[i]
        for v in obj.values():
            fix_jsonld_ids(v, ids)
    if isinstance(obj, MutableSequence):
        for entry in obj:
            fix_jsonld_ids(entry, ids)


def makerdf(
    workflow: Optional[str],
    wf: Union[CommentedMap, float, str, CommentedSeq],
    ctx: ContextType,
    graph: Optional[Graph] = None,
) -> Graph:
    prefixes: Final = {}
    idfields: Final = []
    for k, v in ctx.items():
        if isinstance(v, MutableMapping):
            url = v["@id"]
        else:
            url = v
        if url == "@id":
            idfields.append(k)
        doc_url, frg = urldefrag(url)
        if "/" in frg:
            p = frg.split("/")[0]
            prefixes[p] = f"{doc_url}#{p}/"

    fix_jsonld_ids(wf, idfields)

    g: Final = Graph() if graph is None else graph

    if isinstance(wf, MutableSequence):
        for w in wf:
            w["@context"] = ctx
            g.parse(
                data=json_dumps(w, default=str),
                format="json-ld",
                publicID=str(workflow),
            )
    elif isinstance(wf, MutableMapping):
        wf["@context"] = ctx
        g.parse(data=json_dumps(wf, default=str), format="json-ld", publicID=str(workflow))
    else:
        raise SchemaException(f"{wf} is not a workflow")

    # Bug in json-ld loader causes @id fields to be added to the graph
    for sub, pred, obj in g.triples((None, URIRef("@id"), None)):
        g.remove((sub, pred, obj))

    for k2, v2 in prefixes.items():
        g.namespace_manager.bind(k2, v2)

    return g
