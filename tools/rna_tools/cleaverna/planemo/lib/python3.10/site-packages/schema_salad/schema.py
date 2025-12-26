"""Functions to process Schema Salad schemas."""

import copy
import hashlib
from collections.abc import Mapping, MutableMapping, MutableSequence
from importlib.resources import files
from typing import IO, Any, Final, Optional, Union, cast
from urllib.parse import urlparse

from ruamel.yaml.comments import CommentedMap, CommentedSeq

from schema_salad.utils import (
    CacheType,
    ResolveType,
    add_dictlist,
    aslist,
    convert_to_dict,
    flatten,
    json_dumps,
    yaml_no_ts,
)

from . import _logger, jsonld_context, ref_resolver, validate
from .avro.schema import Names, SchemaParseException, is_subtype, make_avsc_object
from .exceptions import (
    ClassValidationException,
    SchemaSaladException,
    ValidationException,
)
from .ref_resolver import Loader
from .sourceline import SourceLine, add_lc_filename, relname

SALAD_FILES: Final = (
    "metaschema.yml",
    "metaschema_base.yml",
    "salad.md",
    "field_name.yml",
    "import_include.md",
    "link_res.yml",
    "ident_res.yml",
    "vocab_res.yml",
    "vocab_res.yml",
    "field_name_schema.yml",
    "field_name_src.yml",
    "field_name_proc.yml",
    "ident_res_schema.yml",
    "ident_res_src.yml",
    "ident_res_proc.yml",
    "link_res_schema.yml",
    "link_res_src.yml",
    "link_res_proc.yml",
    "vocab_res_schema.yml",
    "vocab_res_src.yml",
    "vocab_res_proc.yml",
    "map_res.yml",
    "map_res_schema.yml",
    "map_res_src.yml",
    "map_res_proc.yml",
    "typedsl_res.yml",
    "typedsl_res_schema.yml",
    "typedsl_res_src.yml",
    "typedsl_res_proc.yml",
    "sfdsl_res.yml",
    "sfdsl_res_schema.yml",
    "sfdsl_res_src.yml",
    "sfdsl_res_proc.yml",
)

saladp: Final = "https://w3id.org/cwl/salad#"
primitives: Final = {
    "http://www.w3.org/2001/XMLSchema#string",
    "http://www.w3.org/2001/XMLSchema#boolean",
    "http://www.w3.org/2001/XMLSchema#int",
    "http://www.w3.org/2001/XMLSchema#long",
    saladp + "null",
    saladp + "enum",
    saladp + "array",
    saladp + "record",
    saladp + "Any",
}


cached_metaschema: Optional[tuple[Names, list[dict[str, str]], Loader]] = None


def get_metaschema() -> tuple[Names, list[dict[str, str]], Loader]:
    """Instantiate the metaschema."""
    global cached_metaschema
    if cached_metaschema is not None:
        return cached_metaschema

    loader: Final = ref_resolver.Loader(
        {
            "Any": saladp + "Any",
            "ArraySchema": saladp + "ArraySchema",
            "Array_symbol": saladp + "ArraySchema/type/Array_symbol",
            "DocType": saladp + "DocType",
            "Documentation": saladp + "Documentation",
            "Documentation_symbol": saladp + "Documentation/type/Documentation_symbol",
            "Documented": saladp + "Documented",
            "EnumSchema": saladp + "EnumSchema",
            "Enum_symbol": saladp + "EnumSchema/type/Enum_symbol",
            "JsonldPredicate": saladp + "JsonldPredicate",
            "MapSchema": saladp + "MapSchema",
            "Map_symbol": saladp + "MapSchema/type/Map_symbol",
            "NamedType": saladp + "NamedType",
            "PrimitiveType": saladp + "PrimitiveType",
            "RecordField": saladp + "RecordField",
            "RecordSchema": saladp + "RecordSchema",
            "Record_symbol": saladp + "RecordSchema/type/Record_symbol",
            "SaladEnumSchema": saladp + "SaladEnumSchema",
            "SaladRecordField": saladp + "SaladRecordField",
            "SaladRecordSchema": saladp + "SaladRecordSchema",
            "SchemaDefinedType": saladp + "SchemaDefinedType",
            "SpecializeDef": saladp + "SpecializeDef",
            "UnionSchema": saladp + "UnionSchema",
            "Union_symbol": saladp + "UnionSchema/type/Union_symbol",
            "_container": saladp + "JsonldPredicate/_container",
            "_id": {"@id": saladp + "_id", "@type": "@id", "identity": True},
            "_type": saladp + "JsonldPredicate/_type",
            "abstract": saladp + "SaladRecordSchema/abstract",
            "array": saladp + "array",
            "boolean": "http://www.w3.org/2001/XMLSchema#boolean",
            "dct": "http://purl.org/dc/terms/",
            "default": {"@id": saladp + "default", "noLinkCheck": True},
            "doc": "rdfs:comment",
            "docAfter": {"@id": saladp + "docAfter", "@type": "@id"},
            "docChild": {"@id": saladp + "docChild", "@type": "@id"},
            "docParent": {"@id": saladp + "docParent", "@type": "@id"},
            "documentRoot": saladp + "SchemaDefinedType/documentRoot",
            "documentation": saladp + "documentation",
            "double": "http://www.w3.org/2001/XMLSchema#double",
            "enum": saladp + "enum",
            "extends": {"@id": saladp + "extends", "@type": "@id", "refScope": 1},
            "fields": {
                "@id": saladp + "fields",
                "mapPredicate": "type",
                "mapSubject": "name",
            },
            "float": "http://www.w3.org/2001/XMLSchema#float",
            "identity": saladp + "JsonldPredicate/identity",
            "inVocab": saladp + "NamedType/inVocab",
            "int": "http://www.w3.org/2001/XMLSchema#int",
            "items": {"@id": saladp + "items", "@type": "@vocab", "refScope": 2},
            "jsonldPredicate": "sld:jsonldPredicate",
            "long": "http://www.w3.org/2001/XMLSchema#long",
            "map": saladp + "map",
            "mapPredicate": saladp + "JsonldPredicate/mapPredicate",
            "mapSubject": saladp + "JsonldPredicate/mapSubject",
            "name": "@id",
            "names": {"@id": saladp + "names", "@type": "@vocab", "refScope": 2},
            "noLinkCheck": saladp + "JsonldPredicate/noLinkCheck",
            "null": saladp + "null",
            "rdf": "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
            "rdfs": "http://www.w3.org/2000/01/rdf-schema#",
            "record": saladp + "record",
            "refScope": saladp + "JsonldPredicate/refScope",
            "sld": saladp,
            "specialize": {
                "@id": saladp + "specialize",
                "mapPredicate": "specializeTo",
                "mapSubject": "specializeFrom",
            },
            "specializeFrom": {
                "@id": saladp + "specializeFrom",
                "@type": "@id",
                "refScope": 1,
            },
            "specializeTo": {
                "@id": saladp + "specializeTo",
                "@type": "@id",
                "refScope": 1,
            },
            "string": "http://www.w3.org/2001/XMLSchema#string",
            "subscope": saladp + "JsonldPredicate/subscope",
            "symbols": {"@id": saladp + "symbols", "@type": "@id", "identity": True},
            "type": {
                "@id": saladp + "type",
                "@type": "@vocab",
                "refScope": 2,
                "typeDSL": True,
            },
            "typeDSL": saladp + "JsonldPredicate/typeDSL",
            "union": saladp + "union",
            "values": {"@id": saladp + "values", "@type": "@vocab", "refScope": 2},
            "xsd": "http://www.w3.org/2001/XMLSchema#",
        },
        salad_version="v1.3",
    )

    for salad in SALAD_FILES:
        loader.cache["https://w3id.org/cwl/" + salad] = (
            files("schema_salad").joinpath("metaschema/" + salad).read_text("UTF-8")
        )

    loader.cache["https://w3id.org/cwl/salad"] = (
        files("schema_salad").joinpath("metaschema/metaschema.yml").read_text("UTF-8")
    )

    yaml: Final = yaml_no_ts()
    j: Final = yaml.load(loader.cache["https://w3id.org/cwl/salad"])
    add_lc_filename(j, "metaschema.yml")
    j2: Final = loader.resolve_all(j, saladp)[0]

    if not isinstance(j2, list):
        _logger.error("%s", j2)
        raise SchemaParseException(f"Not a list: {j2}")
    sch_obj: Final = make_avro(j2, loader, loader.vocab)
    try:
        sch_names: Final = make_avro_schema_from_avro(sch_obj)
    except SchemaParseException:
        _logger.error("Metaschema error, avro was:\n%s", json_dumps(sch_obj, indent=4))
        raise
    validate_doc(sch_names, j2, loader, strict=True)
    cached_metaschema = (sch_names, j2, loader)
    return cached_metaschema


def add_namespaces(metadata: Mapping[str, Any], namespaces: MutableMapping[str, str]) -> None:
    """Collect the provided namespaces, checking for conflicts."""
    for key, value in metadata.items():
        if key not in namespaces:
            namespaces[key] = value
        elif namespaces[key] != value:
            raise ValidationException(
                f"Namespace prefix {key!r} has conflicting definitions {namespaces[key]!r}"
                " and {value!r}."
            )


def collect_namespaces(metadata: Mapping[str, Any]) -> dict[str, str]:
    """Walk through the metadata object, collecting namespace declarations."""
    namespaces: dict[str, str] = {}
    if "$import_metadata" in metadata:
        for value in metadata["$import_metadata"].values():
            add_namespaces(collect_namespaces(value), namespaces)
    if "$namespaces" in metadata:
        add_namespaces(metadata["$namespaces"], namespaces)
    return namespaces


schema_type = tuple[Loader, Union[Names, SchemaParseException], dict[str, Any], Loader]


def load_schema(
    schema_ref: ResolveType,
    cache: Optional[CacheType] = None,
) -> schema_type:
    """
    Load a schema that can be used to validate documents using load_and_validate.

    :returns: document_loader, avsc_names, schema_metadata, metaschema_loader
    """
    metaschema_names, _metaschema_doc, metaschema_loader = get_metaschema()
    if cache is not None:
        # we want to replace some items in the cache, so we need to
        # make a new Loader with an empty index.
        for k, v in metaschema_loader.cache.items():
            if k not in cache:
                cache[k] = v
        metaschema_loader = Loader(
            ctx=metaschema_loader.ctx, cache=cache, session=metaschema_loader.session
        )
    schema_doc, schema_metadata = metaschema_loader.resolve_ref(schema_ref, "")

    if not isinstance(schema_doc, MutableSequence):
        raise ValidationException("Schema reference must resolve to a list.")

    validate_doc(metaschema_names, schema_doc, metaschema_loader, True)
    metactx = schema_metadata.get("@context", {})
    metactx.update(collect_namespaces(schema_metadata))
    schema_ctx = jsonld_context.salad_to_jsonld_context(schema_doc, metactx)[0]

    # Create the loader that will be used to load the target document.
    document_loader = Loader(schema_ctx, cache=cache)

    # Make the Avro validation that will be used to validate the target
    # document
    avsc_names = make_avro_schema(schema_doc, document_loader, metaschema_loader.vocab)

    return document_loader, avsc_names, schema_metadata, metaschema_loader


def load_and_validate(
    document_loader: Loader,
    avsc_names: Names,
    document: Union[CommentedMap, str],
    strict: bool,
    strict_foreign_properties: bool = False,
) -> tuple[Any, dict[str, Any]]:
    """Load a document and validate it with the provided schema.

    return data, metadata
    """
    try:
        if isinstance(document, CommentedMap):
            data, metadata = document_loader.resolve_all(
                document,
                document["id"],
                checklinks=True,
                strict_foreign_properties=strict_foreign_properties,
            )
        else:
            data, metadata = document_loader.resolve_ref(
                document,
                checklinks=True,
                strict_foreign_properties=strict_foreign_properties,
            )

        validate_doc(
            avsc_names,
            data,
            document_loader,
            strict,
            strict_foreign_properties=strict_foreign_properties,
        )
    except ValidationException as exc:
        raise ValidationException("", None, [exc]) from exc
    return data, metadata


def validate_doc(
    schema_names: Names,
    doc: ResolveType,
    loader: Loader,
    strict: bool,
    strict_foreign_properties: bool = False,
) -> None:
    """Validate a document using the provided schema."""
    has_root = False
    for root in schema_names.names.values():
        if (hasattr(root, "get_prop") and root.get_prop("documentRoot")) or (
            "documentRoot" in root.props
        ):
            has_root = True
            break

    if not has_root:
        raise ValidationException("No document roots defined in the schema")

    if isinstance(doc, MutableSequence):
        vdoc = doc
    elif isinstance(doc, CommentedMap):
        vdoc = CommentedSeq([doc])
        vdoc.lc.add_kv_line_col(0, [doc.lc.line, doc.lc.col])
        vdoc.lc.filename = doc.lc.filename
    else:
        raise ValidationException("Document must be dict or list")

    roots: Final = []
    for root in schema_names.names.values():
        if (hasattr(root, "get_prop") and root.get_prop("documentRoot")) or (
            root.props.get("documentRoot")
        ):
            roots.append(root)

    anyerrors: Final = []
    for pos, item in enumerate(vdoc):
        sourceline = SourceLine(vdoc, pos, str)
        success = False
        for root in roots:
            success = validate.validate_ex(
                root,
                item,
                loader.identifiers,
                strict,
                foreign_properties=loader.foreign_properties,
                raise_ex=False,
                skip_foreign_properties=loader.skip_schemas,
                strict_foreign_properties=strict_foreign_properties,
                vocab=loader.vocab,
            )
            if success:
                break

        if not success:
            errors: list[SchemaSaladException] = []
            for root in roots:
                if hasattr(root, "get_prop"):
                    name = root.get_prop("name")
                elif hasattr(root, "name"):
                    name = root.name

                try:
                    validate.validate_ex(
                        root,
                        item,
                        loader.identifiers,
                        strict,
                        foreign_properties=loader.foreign_properties,
                        raise_ex=True,
                        skip_foreign_properties=loader.skip_schemas,
                        strict_foreign_properties=strict_foreign_properties,
                        vocab=loader.vocab,
                    )
                except ClassValidationException as exc1:
                    errors = [
                        ClassValidationException(
                            f"tried {validate.friendly(name)!r} but", sourceline, [exc1]
                        )
                    ]
                    break
                except ValidationException as exc2:
                    errors.append(
                        ValidationException(
                            f"tried {validate.friendly(name)!r} but", sourceline, [exc2]
                        )
                    )

            objerr = "Invalid"
            for ident in loader.identifiers:
                if ident in item:
                    objerr = f"Object {relname(item[ident])!r} is not valid because"
                    break
            anyerrors.append(ValidationException(objerr, sourceline, errors, "-"))
    if anyerrors:
        raise ValidationException("", None, anyerrors, "*")


def get_anon_name(rec: MutableMapping[str, Union[str, dict[str, str], list[str]]]) -> str:
    """Calculate a reproducible name for anonymous types."""
    if "name" in rec:
        name: Final = rec["name"]
        if isinstance(name, str):
            return name
        raise ValidationException(f"Expected name field to be a string, was {name}")
    anon_name = ""
    if rec["type"] in ("enum", saladp + "enum"):
        for sym in rec["symbols"]:
            anon_name += sym
        return "anon.enum_" + hashlib.sha1(anon_name.encode("UTF-8")).hexdigest()  # nosec
    if rec["type"] in ("record", saladp + "record"):
        for field in rec["fields"]:
            if isinstance(field, Mapping):
                anon_name += field["name"]
            else:
                raise ValidationException(
                    f"Expected entries in 'fields' to also be maps, was {field}."
                )
        return "record_" + hashlib.sha1(anon_name.encode("UTF-8")).hexdigest()  # nosec
    if rec["type"] in ("array", saladp + "array"):
        return ""
    raise ValidationException("Expected enum or record, was {rec['type'])}")


def replace_type(
    items: Any,
    spec: dict[str, Any],
    loader: Loader,
    found: set[str],
    find_embeds: bool = True,
    deepen: bool = True,
) -> Any:
    """Go through and replace types in the 'spec' mapping."""
    if isinstance(items, MutableMapping):
        # recursively check these fields for types to replace
        if items.get("type") in ("record", "enum", "map", "union") and items.get("name"):
            if items["name"] in found:
                return items["name"]
            found.add(items["name"])

        if not deepen:
            return items

        items = copy.copy(items)
        if not items.get("name"):
            items["name"] = get_anon_name(items)
        for name in ("type", "items", "fields", "values"):
            if name in items:
                items[name] = replace_type(
                    items[name],
                    spec,
                    loader,
                    found,
                    find_embeds=find_embeds,
                    deepen=find_embeds,
                )
                if isinstance(items[name], MutableSequence):
                    items[name] = flatten(items[name])

        return items
    if isinstance(items, MutableSequence):
        # recursively transform list
        return [
            replace_type(i, spec, loader, found, find_embeds=find_embeds, deepen=deepen)
            for i in items
        ]
    if isinstance(items, str):
        # found a string which is a symbol corresponding to a type.
        replace_with = None
        if items in loader.vocab:
            # If it's a vocabulary term, first expand it to its fully qualified
            # URI
            items = loader.vocab[items]

        if items in spec:
            # Look up in specialization map
            replace_with = spec[items]

        if replace_with:
            return replace_type(replace_with, spec, loader, found, find_embeds=find_embeds)
        found.add(items)
    return items


def avro_field_name(url: str) -> str:
    """
    Turn a URL into an Avro-safe name.

    If the URL has no fragment, return this plain URL.

    Extract either the last part of the URL fragment past the slash, otherwise
    the whole fragment.
    """
    d: Final = urlparse(url)
    if d.fragment:
        return d.fragment.split("/")[-1]
    return d.path.split("/")[-1]


Avro = Union[MutableMapping[str, Any], MutableSequence[Any], str]


def make_valid_avro(
    items: Avro,
    alltypes: dict[str, dict[str, Any]],
    found: set[str],
    union: bool = False,
    fielddef: bool = False,
    vocab: Optional[dict[str, str]] = None,
) -> Union[Avro, MutableMapping[str, str], str, list[Union[Any, MutableMapping[str, str], str]]]:
    """Convert our schema to be more avro like."""
    if vocab is None:
        _, _, metaschema_loader = get_metaschema()
        vocab = metaschema_loader.vocab

    # Possibly could be integrated into our fork of avro/schema.py?
    if isinstance(items, MutableMapping):
        avro: Final = copy.copy(items)
        if avro.get("name"):
            if fielddef:
                avro["name"] = avro_field_name(avro["name"])
            else:
                avro["name"] = validate.avro_type_name(avro["name"])

        if "type" in avro and avro["type"] in (
            saladp + "record",
            saladp + "enum",
            saladp + "map",
            saladp + "union",
            "record",
            "enum",
            "map",
            "union",
        ):
            if (hasattr(avro, "get") and avro.get("abstract")) or ("abstract" in avro):
                return avro
            if "name" in avro:
                if avro["name"] in found:
                    return cast(str, avro["name"])
                found.add(avro["name"])
        for field in ("type", "items", "names", "values", "fields"):
            if field in avro:
                avro[field] = make_valid_avro(
                    avro[field],
                    alltypes,
                    found,
                    union=True,
                    fielddef=(field == "fields"),
                    vocab=vocab,
                )
        if "symbols" in avro:
            avro["symbols"] = [avro_field_name(sym) for sym in avro["symbols"]]
        return avro
    if items and isinstance(items, MutableSequence):
        ret: Final = []
        for i in items:
            ret.append(
                make_valid_avro(i, alltypes, found, union=union, fielddef=fielddef, vocab=vocab)
            )
        return ret
    if union and isinstance(items, str):
        if items in alltypes and validate.avro_type_name(items) not in found:
            return make_valid_avro(alltypes[items], alltypes, found, union=union, vocab=vocab)
        if items in vocab:
            return validate.avro_type_name(vocab[items])
        return validate.avro_type_name(items)
    return items


def deepcopy_strip(item: Any) -> Any:
    """
    Make a deep copy of list and dict objects.

    Intentionally do not copy attributes.  This is to discard CommentedMap and
    CommentedSeq metadata which is very expensive with regular copy.deepcopy.
    """
    if isinstance(item, MutableMapping):
        return {k: deepcopy_strip(v) for k, v in item.items()}
    if isinstance(item, MutableSequence):
        return [deepcopy_strip(k) for k in item]
    return item


def extend_and_specialize(items: list[dict[str, Any]], loader: Loader) -> list[dict[str, Any]]:
    """Apply 'extend' and 'specialize' to fully materialize derived record types."""
    items2: Final = deepcopy_strip(items)
    types: dict[str, Any] = {
        i["name"]: i for i in items2
    }  # no Final, error: ‘CPyStatic_types___’ undeclared (first use in this function)
    types.update({k[len(saladp) :]: v for k, v in types.items() if k.startswith(saladp)})
    results: Final = []

    for stype in items2:
        if "extends" in stype:
            specs: dict[str, str] = {}
            if "specialize" in stype:
                for spec in aslist(stype["specialize"]):
                    specs[spec["specializeFrom"]] = spec["specializeTo"]

            exfields: list[Any] = []
            exsym: list[str] = []
            for ex in aslist(stype["extends"]):
                if ex not in types:
                    raise ValidationException(
                        f"Extends {stype['extends']} in {stype['name']} refers to invalid base type."
                    )

                basetype = copy.copy(types[ex])

                if stype["type"] == "record":
                    if specs:
                        basetype["fields"] = replace_type(
                            basetype.get("fields", []), specs, loader, set()
                        )

                    for field in basetype.get("fields", []):
                        if "inherited_from" not in field:
                            field["inherited_from"] = ex

                    exfields.extend(basetype.get("fields", []))
                elif stype["type"] == "enum":
                    exsym.extend(basetype.get("symbols", []))

            if stype["type"] == "record":
                stype = copy.copy(stype)
                combined_fields = []
                fields = stype.get("fields", [])
                # We use short names here so that if a type inherits a field
                # (e.g. Child#id) from a parent (Parent#id) we avoid adding
                # the same field twice (previously we had just
                # ``exfields.extends(stype.fields)``).
                sns_fields = {shortname(field["name"]): field for field in fields}
                sns_exfields = {shortname(exfield["name"]): exfield for exfield in exfields}

                # N.B.: This could be simpler. We could have a single loop
                #       to create the list of fields. The reason for this more
                #       convoluted solution is to make sure we keep the order
                #       of ``exfields`` first, and then the type fields. That's
                #       because we have unit tests that rely on the order that
                #       fields are written. Codegen output changes as well.
                #       We are relying on the insertion order preserving
                #       property of Python dicts (i.e. relying on Py3.5+).

                # First pass adding the exfields.
                for sn_exfield, exfield in sns_exfields.items():
                    field = sns_fields.get(sn_exfield, None)
                    if field is None:
                        field = exfield
                    else:
                        # make sure field name has not been used yet
                        if not is_subtype(types, exfield["type"], field["type"]):
                            raise SchemaParseException(
                                f"Field name {field['name']} already in use with "
                                "incompatible type. "
                                f"{field['type']} vs {exfield['type']}."
                            )
                    combined_fields.append(field)
                # Second pass, now add the ones that are specific to the subtype.
                for field in sns_fields.values():
                    if field not in combined_fields:
                        combined_fields.append(field)

                stype["fields"] = combined_fields

                fieldnames: set[str] = set()
                for field in stype["fields"]:
                    if field["name"] in fieldnames:
                        raise ValidationException(
                            f"Field name {field['name']} appears twice in {stype['name']}"
                        )
                    fieldnames.add(field["name"])
            elif stype["type"] == "enum":
                stype = copy.copy(stype)
                exsym.extend(stype.get("symbols", []))
                stype["symbols"] = exsym

            types[stype["name"]] = stype

        results.append(stype)

    ex_types: Final = {}
    for result in results:
        ex_types[result["name"]] = result

    extended_by: Final[dict[str, str]] = {}
    for result in results:
        if "extends" in result:
            for ex in aslist(result["extends"]):
                if ex_types[ex].get("abstract"):
                    add_dictlist(extended_by, ex, ex_types[result["name"]])
                    add_dictlist(extended_by, validate.avro_type_name(ex), ex_types[ex])

    for result in results:
        if result.get("abstract") and result["name"] not in extended_by:
            raise ValidationException(
                f"{result['name']} is abstract but missing a concrete subtype"
            )

    for result in results:
        if "fields" in result:
            result["fields"] = replace_type(result["fields"], extended_by, loader, set())
        elif "values" in result:
            result["values"] = replace_type(result["values"], extended_by, loader, set())

    return results


def make_avro(
    i: list[dict[str, Any]],
    loader: Loader,
    metaschema_vocab: Optional[dict[str, str]] = None,
) -> list[Any]:
    j: Final = extend_and_specialize(i, loader)

    name_dict: Final[dict[str, dict[str, Any]]] = {}
    for entry in j:
        name_dict[entry["name"]] = entry

    avro: Final = make_valid_avro(j, name_dict, set(), vocab=metaschema_vocab)

    return [
        t
        for t in avro
        if isinstance(t, MutableMapping)
        and not t.get("abstract")
        and t.get("type") != "org.w3id.cwl.salad.documentation"
    ]


def make_avro_schema(
    i: list[Any], loader: Loader, metaschema_vocab: Optional[dict[str, str]] = None
) -> Names:
    """
    All in one convenience function.

    Call make_avro() and make_avro_schema_from_avro() separately if you need
    the intermediate result for diagnostic output.
    """
    names: Final = Names()
    avro: Final = make_avro(i, loader, metaschema_vocab)
    make_avsc_object(convert_to_dict(avro), names)
    return names


def make_avro_schema_from_avro(avro: list[Union[Avro, dict[str, str], str]]) -> Names:
    """Create avro.schema.Names from the given definitions."""
    names: Final = Names()
    make_avsc_object(convert_to_dict(avro), names)
    return names


def shortname(inputid: str) -> str:
    """Return the last segment of the provided fragment or path."""
    parsed_id: Final = urlparse(inputid)
    if parsed_id.fragment:
        return parsed_id.fragment.split("/")[-1]
    return parsed_id.path.split("/")[-1]


def print_inheritance(doc: list[dict[str, Any]], stream: IO[Any]) -> None:
    """Write a Grapviz inheritance graph for the supplied document."""
    stream.write("digraph {\n")
    for entry in doc:
        if entry["type"] == "record":
            label = name = shortname(entry["name"])
            fields = entry.get("fields", [])
            if fields:
                label += "\\n* {}\\l".format(
                    "\\l* ".join(shortname(field["name"]) for field in fields)
                )
            shape = "ellipse" if entry.get("abstract") else "box"
            stream.write(f'"{name}" [shape={shape} label="{label}"];\n')  # noqa: B907
            if "extends" in entry:
                for target in aslist(entry["extends"]):
                    stream.write(f'"{shortname(target)}" -> "{name}";\n')  # noqa: B907
    stream.write("}\n")


def print_fieldrefs(doc: list[dict[str, Any]], loader: Loader, stream: IO[Any]) -> None:
    """Write a GraphViz graph of the relationships between the fields."""
    obj: Final = extend_and_specialize(doc, loader)

    stream.write("digraph {\n")
    for entry in obj:
        if entry.get("abstract"):
            continue
        if entry["type"] == "record":
            label = shortname(entry["name"])
            for field in entry.get("fields", []):
                found: set[str] = set()
                field_name = shortname(field["name"])
                replace_type(field["type"], {}, loader, found, find_embeds=False)
                for each_type in found:
                    if each_type not in primitives:
                        stream.write(
                            f"{label!r} -> {shortname(each_type)!r} [label={field_name!r}];\n"
                        )
    stream.write("}\n")
