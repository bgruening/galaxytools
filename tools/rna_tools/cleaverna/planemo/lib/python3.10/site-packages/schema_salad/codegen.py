"""Generate language specific loaders for a particular SALAD schema."""

import sys
from collections.abc import MutableMapping, MutableSequence
from io import TextIOWrapper
from typing import Any, Final, Optional, TextIO, Union
from urllib.parse import urlsplit

from . import schema
from .codegen_base import CodeGenBase
from .cpp_codegen import CppCodeGen
from .dlang_codegen import DlangCodeGen
from .dotnet_codegen import DotNetCodeGen
from .exceptions import SchemaSaladException
from .java_codegen import JavaCodeGen
from .python_codegen import PythonCodeGen
from .ref_resolver import Loader
from .schema import shortname
from .typescript_codegen import TypeScriptCodeGen
from .utils import aslist

FIELD_SORT_ORDER: Final = ["class", "id", "name"]


def codegen(
    lang: str,
    i: list[dict[str, str]],
    schema_metadata: dict[str, Any],
    loader: Loader,
    target: Optional[str] = None,
    examples: Optional[str] = None,
    package: Optional[str] = None,
    copyright: Optional[str] = None,
    spdx_copyright_text: Optional[list[str]] = None,
    spdx_license_identifier: Optional[str] = None,
    parser_info: Optional[str] = None,
) -> None:
    """Generate classes with loaders for the given Schema Salad description."""
    j = schema.extend_and_specialize(i, loader)

    gen: Optional[CodeGenBase] = None
    base = schema_metadata.get("$base", schema_metadata.get("id"))
    # ``urlsplit`` decides whether to return an encoded result based
    # on the object type. To ensure the code behaves the same for Py
    # 3.6+, we enforce that the input value is of type ``str``.
    if base is None:
        base = ""
    sp = urlsplit(base)
    pkg = (
        package
        if package
        else ".".join(list(reversed(sp.netloc.split("."))) + sp.path.strip("/").split("/"))
    )
    info = parser_info or pkg
    salad_version = schema_metadata.get("saladVersion", "v1.1")

    if lang in {"python", "cpp", "dlang"}:
        if target:
            dest: Union[TextIOWrapper, TextIO] = open(target, mode="w", encoding="utf-8")
        else:
            dest = sys.stdout
        if lang == "cpp":
            gen = CppCodeGen(
                base,
                dest,
                examples,
                pkg,
                copyright,
                spdx_copyright_text,
                spdx_license_identifier,
            )
            gen.parse(j)
            return
        if lang == "dlang":
            gen = DlangCodeGen(
                base,
                dest,
                examples,
                pkg,
                copyright,
                info,
                salad_version,
            )
            gen.parse(j)
            return
        gen = PythonCodeGen(
            dest, copyright=copyright, parser_info=info, salad_version=salad_version
        )

    elif lang == "java":
        gen = JavaCodeGen(
            base,
            target=target,
            examples=examples,
            package=pkg,
            copyright=copyright,
        )
    elif lang == "typescript":
        gen = TypeScriptCodeGen(base, target=target, package=pkg, examples=examples)
    elif lang == "dotnet":
        gen = DotNetCodeGen(base, target=target, package=pkg, examples=examples)
    else:
        raise SchemaSaladException(f"Unsupported code generation language {lang!r}")

    gen.prologue()

    document_roots = []

    for rec in j:
        if rec["type"] in ("enum", "map", "record", "union"):
            jld = rec.get("jsonldPredicate")
            if isinstance(jld, MutableMapping):
                gen.type_loader(rec, jld.get("_container"), jld.get("noLinkCheck"))
            else:
                gen.type_loader(rec)
            gen.add_vocab(shortname(rec["name"]), rec["name"])

    for rec in j:
        if rec["type"] == "enum":
            for symbol in rec["symbols"]:
                gen.add_vocab(shortname(symbol), symbol)

        if rec["type"] == "record":
            if rec.get("documentRoot"):
                document_roots.append(rec["name"])

            field_names = []
            optional_fields = set()
            for field in rec.get("fields", []):
                field_name = shortname(field["name"])
                field_names.append(field_name)
                tp = field["type"]
                if isinstance(tp, MutableSequence) and tp[0] == "https://w3id.org/cwl/salad#null":
                    optional_fields.add(field_name)

            idfield = ""
            for field in rec.get("fields", []):
                if field.get("jsonldPredicate") == "@id":
                    idfield = field.get("name")

            gen.begin_class(
                rec["name"],
                aslist(rec.get("extends", [])),
                rec.get("doc", ""),
                rec.get("abstract", False),
                field_names,
                idfield,
                optional_fields,
            )
            gen.add_vocab(shortname(rec["name"]), rec["name"])

            sorted_fields = sorted(
                rec.get("fields", []),
                key=lambda i: (
                    FIELD_SORT_ORDER.index(i["name"].split("/")[-1])
                    if i["name"].split("/")[-1] in FIELD_SORT_ORDER
                    else 100
                ),
            )

            for field in sorted_fields:
                if field.get("jsonldPredicate") == "@id":
                    subscope = field.get("subscope")
                    fieldpred = field["name"]
                    optional = bool("https://w3id.org/cwl/salad#null" in field["type"])
                    uri_loader = gen.uri_loader(gen.type_loader(field["type"]), True, False, None)
                    gen.declare_id_field(
                        fieldpred,
                        uri_loader,
                        field.get("doc"),
                        optional,
                    )
                    break

            for field in sorted_fields:
                optional = bool("https://w3id.org/cwl/salad#null" in field["type"])
                jld = field.get("jsonldPredicate")
                fieldpred = field["name"]
                subscope = None

                if isinstance(jld, MutableMapping):
                    type_loader = gen.type_loader(
                        field["type"], jld.get("_container"), jld.get("noLinkCheck")
                    )
                    ref_scope = jld.get("refScope")
                    subscope = jld.get("subscope")
                    if jld.get("typeDSL"):
                        type_loader = gen.typedsl_loader(type_loader, ref_scope)
                    elif jld.get("secondaryFilesDSL"):
                        type_loader = gen.secondaryfilesdsl_loader(type_loader)
                    elif jld.get("_type") == "@id":
                        type_loader = gen.uri_loader(
                            type_loader,
                            jld.get("identity", False),
                            False,
                            ref_scope,
                            jld.get("noLinkCheck"),
                        )
                    elif jld.get("_type") == "@vocab":
                        type_loader = gen.uri_loader(
                            type_loader, False, True, ref_scope, jld.get("noLinkCheck")
                        )

                    map_subject = jld.get("mapSubject")
                    if map_subject:
                        type_loader = gen.idmap_loader(
                            field["name"],
                            type_loader,
                            map_subject,
                            jld.get("mapPredicate"),
                        )

                    if "_id" in jld and jld["_id"][0] != "@":
                        fieldpred = jld["_id"]
                else:
                    type_loader = gen.type_loader(field["type"])

                if jld == "@id":
                    continue

                gen.declare_field(fieldpred, type_loader, field.get("doc"), optional, subscope)

            gen.end_class(rec["name"], field_names)

    root_type = list(document_roots)
    root_type.append({"type": "array", "items": document_roots})

    gen.epilogue(gen.type_loader(root_type))
