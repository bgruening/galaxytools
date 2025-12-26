"""TypeScript code generator for a given schema salad definition."""

import os
import shutil
import string
from collections.abc import MutableMapping, MutableSequence
from importlib.resources import files
from io import StringIO
from pathlib import Path
from typing import Any, Final, Optional, Union

from . import _logger, schema
from .codegen_base import CodeGenBase, LazyInitDef, TypeDef
from .exceptions import SchemaException
from .java_codegen import _ensure_directory_and_write, _safe_makedirs
from .schema import shortname
from .utils import Traversable


def doc_to_doc_string(doc: Optional[str], indent_level: int = 0) -> str:
    """Generate a documentation string from a schema salad doc field."""
    lead: Final = " " + "  " * indent_level + "* "
    if doc:
        doc_str = "\n".join([f"{lead}{line}" for line in doc.split("\n")])
    else:
        doc_str = ""
    return doc_str


_string_type_def: Final = TypeDef(
    name="strtype",
    init="new _PrimitiveLoader(TypeGuards.String)",
    instance_type="string",
)

_int_type_def: Final = TypeDef(
    name="inttype", init="new _PrimitiveLoader(TypeGuards.Int)", instance_type="number"
)

_float_type_def: Final = TypeDef(
    name="floattype",
    init="new _PrimitiveLoader(TypeGuards.Float)",
    instance_type="number",
)
_bool_type_def: Final = TypeDef(
    name="booltype",
    init="new _PrimitiveLoader(TypeGuards.Bool)",
    instance_type="boolean",
)

_null_type_def: Final = TypeDef(
    name="undefinedtype",
    init="new _PrimitiveLoader(TypeGuards.Undefined)",
    instance_type="undefined",
)

_any_type_def = TypeDef(name="anyType", init="new _AnyLoader()", instance_type="any")

prims: Final = {
    "http://www.w3.org/2001/XMLSchema#string": _string_type_def,
    "http://www.w3.org/2001/XMLSchema#int": _int_type_def,
    "http://www.w3.org/2001/XMLSchema#long": _int_type_def,
    "http://www.w3.org/2001/XMLSchema#float": _float_type_def,
    "http://www.w3.org/2001/XMLSchema#double": _float_type_def,
    "http://www.w3.org/2001/XMLSchema#boolean": _bool_type_def,
    "https://w3id.org/cwl/salad#null": _null_type_def,
    "https://w3id.org/cwl/salad#Any": _any_type_def,
    "string": _string_type_def,
    "int": _int_type_def,
    "long": _int_type_def,
    "float": _float_type_def,
    "double": _float_type_def,
    "boolean": _bool_type_def,
    "null": _null_type_def,
    "Any": _any_type_def,
}


class TypeScriptCodeGen(CodeGenBase):
    """Generation of TypeScript code for a given Schema Salad definition."""

    def __init__(
        self, base: str, examples: Optional[str], target: Optional[str], package: str
    ) -> None:
        """Initialize the TypeScript codegen."""
        super().__init__()
        self.target_dir: Final = Path(target or ".").resolve()
        self.main_src_dir: Final = self.target_dir / "src"
        self.test_resources_dir: Final = self.target_dir / "src" / "test" / "data"
        self.package: Final = package
        self.base_uri: Final = base
        self.record_types: Final[set[str]] = set()
        self.modules: Final[set[str]] = set()
        self.id_field = ""
        self.examples: Final = examples

    def prologue(self) -> None:
        """Trigger to generate the prolouge code."""
        for src_dir in [self.main_src_dir]:
            _safe_makedirs(src_dir)

        for primitive in prims.values():
            self.declare_type(primitive)

    @staticmethod
    def safe_name(name: str) -> str:
        """Generate a safe version of the given name."""
        avn = schema.avro_field_name(name)
        if avn.startswith("anon."):
            avn = avn[5:]
        if avn in (
            "class",
            "in",
            "extends",
            "abstract",
            "default",
            "package",
            "arguments",
        ):
            # reserved words
            avn = avn + "_"

        return avn

    def begin_class(
        self,  # pylint: disable=too-many-arguments
        classname: str,
        extends: MutableSequence[str],
        doc: str,
        abstract: bool,
        field_names: MutableSequence[str],
        idfield: str,
        optional_fields: set[str],
    ) -> None:
        """Produce the header for the given class."""
        self.current_interface = self.safe_name(classname) + "Properties"
        cls: Final = self.safe_name(classname)
        self.current_class = cls
        self.current_class_is_abstract = abstract
        interface_module_name = self.current_interface
        self.current_interface_target_file = self.main_src_dir / f"{interface_module_name}.ts"
        class_module_name = self.current_class
        self.current_class_target_file = self.main_src_dir / f"{class_module_name}.ts"
        self.current_constructor_signature = StringIO()
        self.current_constructor_body = StringIO()
        self.current_loader = StringIO()
        self.current_serializer = StringIO()
        self.current_fieldtypes: dict[str, TypeDef] = {}
        self.idfield = idfield

        doc_string = f"""
/**
 * Auto-generated interface for {classname}
"""
        if doc:
            doc_string += " *\n"
            doc_string += doc_to_doc_string(doc)
            doc_string += "\n"
        doc_string += " */"

        self.record_types.add(f"{self.current_interface}")
        self.modules.add(interface_module_name)
        with open(self.current_interface_target_file, "w") as f:
            _logger.info("Writing file: %s", self.current_interface_target_file)
            if extends:
                ext = "extends Internal." + ", Internal.".join(
                    self.safe_name(e) + "Properties" for e in extends
                )
            else:
                ext = ""
            f.write(
                """
import * as Internal from './util/Internal'

{docstring}
export interface {cls} {ext} {{
                    """.format(
                    docstring=doc_string,
                    cls=f"{self.current_interface}",
                    ext=ext,
                )
            )
        if self.current_class_is_abstract:
            return

        self.record_types.add(cls)
        self.modules.add(class_module_name)
        with open(self.current_interface_target_file, "a") as f:
            f.write(
                """
  extensionFields?: Internal.Dictionary<any>
"""
            )
        doc_string = f"""
/**
 * Auto-generated class implementation for {classname}
"""
        if doc:
            doc_string += " *\n"
            doc_string += doc_to_doc_string(doc)
            doc_string += "\n"
        doc_string += " */"
        with open(self.current_class_target_file, "w") as f:
            _logger.info("Writing file: %s", self.current_class_target_file)
            f.write(
                """
import {{
  Dictionary,
  expandUrl,
  loadField,
  LoaderInstances,
  LoadingOptions,
  Saveable,
  ValidationException,
  prefixUrl,
  save,
  saveRelativeUri
}} from './util/Internal'
import {{ v4 as uuidv4 }} from 'uuid'
import * as Internal from './util/Internal'

{docstring}
export class {cls} extends Saveable implements Internal.{current_interface} {{
  extensionFields?: Internal.Dictionary<any>
""".format(
                    cls=cls,
                    current_interface=self.current_interface,
                    docstring=doc_string,
                )
            )
        self.current_constructor_signature.write(
            "\n" + "\n" + "  constructor ({loadingOptions, extensionFields"
        )
        self.current_constructor_body.write(
            """
    super(loadingOptions)
    this.extensionFields = extensionFields ?? {}
"""
        )
        self.current_loader.write(
            """
  /**
   * Used to construct instances of {{@link {cls} }}.
   *
   * @param __doc                           Document fragment to load this record object from.
   * @param baseuri                         Base URI to generate child document IDs against.
   * @param loadingOptions                  Context for loading URIs and populating objects.
   * @param docRoot                         ID at this position in the document (if available)
   * @returns                               An instance of {{@link {cls} }}
   * @throws {{@link ValidationException}}    If the document fragment is not a
   *                                        {{@link Dictionary}} or validation of fields fails.
   */
  static override async fromDoc (__doc: any, baseuri: string, loadingOptions: LoadingOptions,
    docRoot?: string): Promise<Saveable> {{
    const _doc = Object.assign({{}}, __doc)
    const __errors: ValidationException[] = []
            """.format(
                cls=cls
            )
        )

        self.current_serializer.write(
            """
  save (top: boolean = false, baseUrl: string = '', relativeUris: boolean = true)
  : Dictionary<any> {
    const r: Dictionary<any> = {}
    for (const ef in this.extensionFields) {
      r[prefixUrl(ef, this.loadingOptions.vocab)] = this.extensionFields.ef
    }
"""
        )

    def end_class(self, classname: str, field_names: list[str]) -> None:
        """Signal that we are done with this class."""
        with open(self.current_interface_target_file, "a") as f:
            f.write("}")
        if self.current_class_is_abstract:
            return

        self.current_constructor_signature.write(
            f"}} : {{loadingOptions?: LoadingOptions}} & Internal.{self.current_interface}) {{"
        )
        self.current_constructor_body.write("  }\n")
        self.current_loader.write(
            """
    const extensionFields: Dictionary<any> = {{}}
    for (const [key, value] of Object.entries(_doc)) {{
      if (!{classname}.attr.has(key)) {{
        if ((key as string).includes(':')) {{
          const ex = expandUrl(key, '', loadingOptions, false, false)
          extensionFields[ex] = value
        }} else {{
          __errors.push(
            new ValidationException(`invalid field ${{key as string}}, \\
            expected one of: {fields}`)
          )
          break
        }}
      }}
    }}

    if (__errors.length > 0) {{
      throw new ValidationException("Trying '{classname}'", __errors)
    }}

    const schema = new {classname}({{
      extensionFields: extensionFields,
      loadingOptions: loadingOptions,
      """.format(
                classname=self.current_class,
                fields=",".join(["\\`" + f + "\\`" for f in field_names]),
            )
        )
        self.current_loader.write(
            ",\n      ".join(self.safe_name(f) + ": " + self.safe_name(f) for f in field_names)
            + "\n    })"
        )
        self.current_loader.write(
            """
    return schema
  }
        """
        )
        self.current_serializer.write(
            """
    if (top) {
      if (this.loadingOptions.namespaces != null) {
        r.$namespaces = this.loadingOptions.namespaces
      }
      if (this.loadingOptions.schemas != null) {
        r.$schemas = this.loadingOptions.schemas
      }
    }
    return r
  }
            """
        )
        with open(
            self.current_class_target_file,
            "a",
        ) as f:
            f.write(self.current_constructor_signature.getvalue())
            f.write(self.current_constructor_body.getvalue())
            f.write(self.current_loader.getvalue())
            f.write(self.current_serializer.getvalue())
            f.write(
                "\n"
                + "  static attr: Set<string> = new Set(["
                + ",".join(["'" + shortname(f) + "'" for f in field_names])
                + "])"
            )
            f.write(
                """
}
"""
            )

    def type_loader(
        self,
        type_declaration: Union[list[Any], dict[str, Any], str],
        container: Optional[str] = None,
        no_link_check: Optional[bool] = None,
    ) -> TypeDef:
        """Parse the given type declaration and declare its components."""
        if isinstance(type_declaration, MutableSequence):
            sub_types1: Final = [self.type_loader(i) for i in type_declaration]
            sub_names: Final[list[str]] = list(dict.fromkeys([i.name for i in sub_types1]))
            sub_instance_types: Final[list[str]] = list(
                dict.fromkeys([i.instance_type for i in sub_types1 if i.instance_type is not None])
            )
            return self.declare_type(
                TypeDef(
                    "unionOf{}".format("Or".join(sub_names)),
                    "new _UnionLoader([{}])".format(", ".join(sub_names)),
                    instance_type=" | ".join(sub_instance_types),
                )
            )
        if isinstance(type_declaration, MutableMapping):
            if type_declaration["type"] in (
                "array",
                "https://w3id.org/cwl/salad#array",
            ):
                i1: Final = self.type_loader(type_declaration["items"])
                return self.declare_type(
                    TypeDef(
                        f"arrayOf{i1.name}",
                        f"new _ArrayLoader([{i1.name}])",
                        instance_type=f"Array<{i1.instance_type}>",
                    )
                )
            if type_declaration["type"] in (
                "map",
                "https://w3id.org/cwl/salad#map",
            ):
                i2: Final = self.type_loader(type_declaration["values"])
                return self.declare_type(
                    TypeDef(
                        f"mapOf{i2.name}",
                        "new _MapLoader([{}], {}, {})".format(
                            i2.name,
                            (
                                f"'{container}'"
                                if container is not None
                                else self.to_typescript(None)
                            ),  # noqa: B907
                            self.to_typescript(no_link_check),
                        ),
                        instance_type=f"Dictionary<{i2.instance_type}>",
                    )
                )
            if type_declaration["type"] in ("enum", "https://w3id.org/cwl/salad#enum"):
                return self.type_loader_enum(type_declaration)

            if type_declaration["type"] in (
                "record",
                "https://w3id.org/cwl/salad#record",
            ):
                return self.declare_type(
                    TypeDef(
                        self.safe_name(type_declaration["name"]) + "Loader",
                        "new _RecordLoader({}.fromDoc, {}, {})".format(
                            self.safe_name(type_declaration["name"]),
                            (
                                f"'{container}'"
                                if container is not None
                                else self.to_typescript(None)
                            ),  # noqa: B907
                            self.to_typescript(no_link_check),
                        ),
                        instance_type="Internal." + self.safe_name(type_declaration["name"]),
                        abstract=type_declaration.get("abstract", False),
                    )
                )
            if type_declaration["type"] in (
                "union",
                "https://w3id.org/cwl/salad#union",
            ):
                # Declare the named loader to handle recursive union definitions
                loader_name: Final = self.safe_name(type_declaration["name"]) + "Loader"
                loader_type: Final = TypeDef(
                    loader_name,
                    "new _UnionLoader([])",
                    instance_type="any",
                )
                self.declare_type(loader_type)
                # Parse inner types
                sub_types2: Final = [self.type_loader(i) for i in type_declaration["names"]]

                # Register lazy initialization for the loader
                self.add_lazy_init(
                    LazyInitDef(
                        loader_name,
                        "{}.addLoaders([{}]);".format(
                            loader_name, ", ".join(s.name for s in sub_types2)
                        ),
                    )
                )
                return loader_type
            raise SchemaException("wft {}".format(type_declaration["type"]))

        if type_declaration in prims:
            return prims[type_declaration]

        if type_declaration in ("Expression", "https://w3id.org/cwl/cwl#Expression"):
            return self.declare_type(
                TypeDef(
                    self.safe_name(type_declaration) + "Loader",
                    "new _ExpressionLoader()",
                    instance_type="string",
                )
            )
        return self.collected_types[self.safe_name(type_declaration) + "Loader"]

    def type_loader_enum(self, type_declaration: dict[str, Any]) -> TypeDef:
        """Build an enum type loader for the given declaration."""
        for sym in type_declaration["symbols"]:
            self.add_vocab(shortname(sym), sym)
        enum_name: Final = self.safe_name(type_declaration["name"])
        enum_module_name: Final = enum_name
        enum_path: Final = self.main_src_dir / f"{enum_module_name}.ts"
        self.modules.add(enum_module_name)
        self.record_types.add(enum_name)
        with open(enum_path, "w") as f:
            _logger.info("Writing file: %s", enum_path)
            f.write(
                """
export enum {enum_name} {{
""".format(
                    enum_name=enum_name
                )
            )
            for sym in type_declaration["symbols"]:
                val = self.safe_name(sym)
                const = self.safe_name(sym).replace("-", "_").replace(".", "_").upper()
                f.write(f"""  {const}='{val}',\n""")  # noqa: B907
            f.write(
                """}
"""
            )
        return self.declare_type(
            TypeDef(
                instance_type="Internal." + enum_name,
                name=self.safe_name(type_declaration["name"]) + "Loader",
                init=f"new _EnumLoader((Object.keys({enum_name}) as Array<keyof typeof "
                f"{enum_name}>).map(key => {enum_name}[key]))",
            )
        )

    def declare_field(
        self,
        name: str,
        fieldtype: TypeDef,
        doc: Optional[str],
        optional: bool,
        subscope: Optional[str],
    ) -> None:
        """Output the code to load the given field."""
        safename: Final = self.safe_name(name)
        fieldname: Final = shortname(name)
        self.current_fieldtypes[safename] = fieldtype
        if fieldtype.instance_type is not None and "undefined" in fieldtype.instance_type:
            optionalstring = "?"
        else:
            optionalstring = ""

        with open(self.current_interface_target_file, "a") as f:
            if doc:
                f.write(
                    """
  /**
{doc_str}
   */
""".format(
                        doc_str=doc_to_doc_string(doc, indent_level=1)
                    )
                )
            if fieldname == "class":
                f.write(
                    "  {safename}{optionalstring}: {type}\n".format(
                        safename=safename,
                        type=fieldtype.instance_type,
                        optionalstring="?",
                    )
                )
            else:
                f.write(
                    "  {safename}{optionalstring}: {type}\n".format(
                        safename=safename,
                        type=fieldtype.instance_type,
                        optionalstring=optionalstring,
                    )
                )
        if self.current_class_is_abstract:
            return

        with open(self.current_class_target_file, "a") as f:
            if doc:
                f.write(
                    """
  /**
{doc_str}
   */
""".format(
                        doc_str=doc_to_doc_string(doc, indent_level=1)
                    )
                )
            f.write(
                "  {safename}{optionalstring}: {type}\n".format(
                    safename=safename,
                    type=fieldtype.instance_type,
                    optionalstring=optionalstring,
                )
            )
        if fieldname == "class":
            if fieldtype.instance_type == "string":
                self.current_constructor_signature.write(
                    f", {safename} = '{self.current_class}'"  # noqa: B907
                )
            else:
                self.current_constructor_signature.write(
                    ", {safename} = {type}.{val}".format(
                        safename=safename,
                        type=fieldtype.instance_type,
                        val=self.current_class.replace("-", "_").replace(".", "_").upper(),
                    )
                )
        else:
            self.current_constructor_signature.write(
                ", {safename}".format(
                    safename=safename,
                )
            )
        self.current_constructor_body.write(
            "    this.{safeName} = {safeName}\n".format(safeName=safename)
        )

        self.current_loader.write(
            """
    let {safename}""".format(
                safename=safename
            )
        )
        if optional:
            self.current_loader.write(
                """
    if ('{fieldname}' in _doc) {{""".format(
                    fieldname=fieldname
                )
            )
            spc = "  "
        else:
            spc = ""

        self.current_loader.write(
            """
{spc}    try {{
{spc}      {safename} = await loadField(_doc.{fieldname}, LoaderInstances.{fieldtype},
{spc}        baseuri, loadingOptions)
{spc}    }} catch (e) {{
{spc}      if (e instanceof ValidationException) {{
{spc}        __errors.push(
{spc}          new ValidationException('the `{fieldname}` field is not valid because: ', [e])
{spc}        )
{spc}      }} else {{
{spc}        throw e
{spc}      }}
{spc}    }}
""".format(
                safename=safename,
                fieldname=fieldname,
                fieldtype=fieldtype.name,
                spc=spc,
            )
        )
        if optional:
            self.current_loader.write("    }\n")

        if name == self.idfield or not self.idfield:
            baseurl = "baseUrl"
        else:
            baseurl = f"this.{self.safe_name(self.idfield)}"

        if fieldtype.is_uri:
            self.current_serializer.write(
                """
    if (this.{safename} != null) {{
      const u = saveRelativeUri(this.{safename}, {base_url}, {scoped_id},
                                relativeUris, {ref_scope})
      if (u != null) {{
        r.{fieldname} = u
      }}
    }}
                """.format(
                    safename=self.safe_name(name),
                    fieldname=shortname(name).strip(),
                    base_url=baseurl,
                    scoped_id=self.to_typescript(fieldtype.scoped_id),
                    ref_scope=self.to_typescript(fieldtype.ref_scope),
                )
            )
        else:
            self.current_serializer.write(
                """
    if (this.{safename} != null) {{
      r.{fieldname} = save(this.{safename}, false, {base_url}, relativeUris)
    }}
                """.format(
                    safename=self.safe_name(name),
                    fieldname=shortname(name).strip(),
                    base_url=baseurl,
                )
            )

    def declare_id_field(
        self,
        name: str,
        fieldtype: TypeDef,
        doc: Optional[str],
        optional: bool,
    ) -> None:
        """Output the code to handle the given ID field."""
        self.declare_field(name, fieldtype, doc, True, "")
        if optional:
            opt = f"""{self.safe_name(name)} = "_" + uuidv4()"""
        else:
            opt = """throw new ValidationException("Missing {fieldname}")""".format(
                fieldname=shortname(name)
            )

        self.current_loader.write(
            """
    const original{safename}IsUndefined = ({safename} === undefined)
    if (original{safename}IsUndefined ) {{
      if (docRoot != null) {{
        {safename} = docRoot
      }} else {{
        {opt}
      }}
    }} else {{
      baseuri = {safename} as string
    }}
            """.format(
                safename=self.safe_name(name), opt=opt
            )
        )

    def to_typescript(self, val: Any) -> Any:
        """Convert a Python keyword to a TypeScript keyword."""
        if val is True:
            return "true"
        elif val is None:
            return "undefined"
        elif val is False:
            return "false"
        return val

    def uri_loader(
        self,
        inner: TypeDef,
        scoped_id: bool,
        vocab_term: bool,
        ref_scope: Optional[int],
        no_link_check: Optional[bool] = None,
    ) -> TypeDef:
        """Construct the TypeDef for the given URI loader."""
        instance_type = inner.instance_type or "any"
        return self.declare_type(
            TypeDef(
                f"uri{inner.name}{scoped_id}{vocab_term}{ref_scope}{no_link_check}",
                "new _URILoader({}, {}, {}, {}, {})".format(
                    inner.name,
                    self.to_typescript(scoped_id),
                    self.to_typescript(vocab_term),
                    self.to_typescript(ref_scope),
                    self.to_typescript(no_link_check),
                ),
                is_uri=True,
                scoped_id=scoped_id,
                ref_scope=ref_scope,
                instance_type=instance_type,
            )
        )

    def idmap_loader(
        self, field: str, inner: TypeDef, map_subject: str, map_predicate: Optional[str]
    ) -> TypeDef:
        """Construct the TypeDef for the given mapped ID loader."""
        instance_type = inner.instance_type or "any"
        return self.declare_type(
            TypeDef(
                f"idmap{self.safe_name(field)}{inner.name}",
                f"new _IdMapLoader({inner.name}, '{map_subject}', '{map_predicate}')",  # noqa: B907
                instance_type=instance_type,
            )
        )

    def typedsl_loader(self, inner: TypeDef, ref_scope: Optional[int]) -> TypeDef:
        """Construct the TypeDef for the given DSL loader."""
        instance_type = inner.instance_type or "any"
        return self.declare_type(
            TypeDef(
                f"typedsl{self.safe_name(inner.name)}{ref_scope}",
                f"new _TypeDSLLoader({self.safe_name(inner.name)}, {ref_scope})",
                instance_type=instance_type,
            )
        )

    def epilogue(self, root_loader: TypeDef) -> None:
        """Trigger to generate the epilouge code."""
        pd = "This project contains TypeScript objects and utilities "
        pd = pd + ' auto-generated by <a href=\\"https://github.com/'
        pd = pd + 'common-workflow-language/schema_salad\\">Schema Salad</a>'
        pd = pd + " for parsing documents corresponding to the "
        pd = pd + str(self.base_uri) + " schema."

        sorted_record_types = sorted(self.record_types)
        generated_class_imports = ",\n  ".join(sorted_record_types)
        template_vars: MutableMapping[str, str] = dict(
            project_name=self.package,
            version="0.0.1-SNAPSHOT",
            project_description=pd,
            license_name="Apache License, Version 2.0",
            generated_class_imports=generated_class_imports,
        )

        def template_from_resource(resource: Traversable) -> string.Template:
            template_str = resource.read_text("utf-8")
            template = string.Template(template_str)
            return template

        def expand_resource_template_to(resource: str, path: Path) -> None:
            template = template_from_resource(
                files("schema_salad").joinpath(f"typescript/{resource}")
            )
            src = template.safe_substitute(template_vars)
            _ensure_directory_and_write(path, src)

        expand_resource_template_to("package.json", self.target_dir / "package.json")
        expand_resource_template_to(".gitignore", self.target_dir / ".gitignore")
        expand_resource_template_to("LICENSE", self.target_dir / "LICENSE")
        expand_resource_template_to("tsconfig.json", self.target_dir / "tsconfig.json")
        expand_resource_template_to("index.ts", self.main_src_dir / "index.ts")

        vocab = ",\n  ".join(
            f"""'{k}': '{self.vocab[k]}'""" for k in sorted(self.vocab.keys())  # noqa: B907
        )
        rvocab = ",\n  ".join(
            f"""'{self.vocab[k]}': '{k}'""" for k in sorted(self.vocab.keys())  # noqa: B907
        )

        loader_instances = ""
        for _, collected_type in self.collected_types.items():
            if not collected_type.abstract:
                loader_instances += "export const {} = {};\n".format(
                    collected_type.name, collected_type.init
                )

        sorted_modules = sorted(self.modules)
        internal_module_exports = "\n".join(f"export * from '../{f}'" for f in sorted_modules)

        if self.lazy_inits:
            loader_instances += "\n"
            for lazy_init in self.lazy_inits.values():
                loader_instances += f"{lazy_init.init}\n"

        example_tests = ""
        if self.examples:
            _safe_makedirs(self.test_resources_dir)
            utils_resources = self.test_resources_dir / "examples"
            if os.path.exists(utils_resources):
                shutil.rmtree(utils_resources)
            shutil.copytree(self.examples, utils_resources)
            for example_name in os.listdir(self.examples):
                if example_name.startswith("valid"):
                    basename = os.path.basename(example_name).rsplit(".", 1)[0]
                    example_tests += """
    it('{basename}', async () => {{
        await loadDocument(__dirname + '/data/examples/{example_name}')
    }})
    it('{basename} by string', async () => {{
        let doc = fs.readFileSync(__dirname + '/data/examples/{example_name}').toString()
        await loadDocumentByString(doc, url.pathToFileURL(__dirname +
            '/data/examples/').toString())
    }})""".format(
                        basename=basename.replace("-", "_").replace(".", "_"),
                        example_name=example_name,
                    )

        template_args: MutableMapping[str, str] = dict(
            internal_module_exports=internal_module_exports,
            loader_instances=loader_instances,
            generated_class_imports=generated_class_imports,
            vocab=vocab,
            rvocab=rvocab,
            root_loader=root_loader.name,
            root_loader_type=root_loader.instance_type or "any",
            tests=example_tests,
        )

        util_src_dirs = {
            "util": self.main_src_dir / "util",
            "test": self.main_src_dir / "test",
        }

        def copy_utils_recursive(util_src: str, util_target: Path) -> None:
            for util in files("schema_salad").joinpath(f"typescript/{util_src}").iterdir():
                if util.is_dir():
                    copy_utils_recursive(os.path.join(util_src, util.name), util_target / util.name)
                    continue
                src_path = util_target / util.name
                src_template = template_from_resource(util)
                src = src_template.safe_substitute(template_args)
                _ensure_directory_and_write(src_path, src)

        for util_src, util_target in util_src_dirs.items():
            copy_utils_recursive(util_src, util_target)

    def secondaryfilesdsl_loader(self, inner: TypeDef) -> TypeDef:
        """Construct the TypeDef for secondary files."""
        instance_type = inner.instance_type or "any"
        return self.declare_type(
            TypeDef(
                f"secondaryfilesdsl{inner.name}",
                f"new _SecondaryDSLLoader({inner.name})",
                instance_type=instance_type,
            )
        )
