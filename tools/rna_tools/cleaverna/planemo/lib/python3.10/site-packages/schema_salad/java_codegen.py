"""Java code generator for a given schema salad definition."""

import os
import re
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
from .schema import shortname
from .utils import Traversable

# experiment at providing more typed objects building a optional type that allows
# referencing one or a list of objects. It is useful for improving the RootLoader
# for simple schema with a single root loader - but doesn't help with CWL at all and
# may even confuse things a bit so turning these off be default.
USE_ONE_OR_LIST_OF_TYPES: Final = False

BASIC_JAVA_IDENTIFIER_RE: Final = re.compile(r"[^0-9a-zA-Z]+")


def _ensure_directory_and_write(path: Path, contents: str) -> None:
    _safe_makedirs(path.parent)
    with open(path, mode="w", encoding="utf-8") as f:
        _logger.info("Writing file: %s", path)
        f.write(contents)


def _doc_to_doc_string(doc: Optional[str], indent_level: int = 0) -> str:
    lead: Final = " " + "  " * indent_level + "* " * indent_level
    if doc:
        doc_str = f"{lead}<BLOCKQUOTE>\n"
        doc_str += "\n".join([f"{lead}{line}" for line in doc.split("\n")])
        doc_str += f"{lead}</BLOCKQUOTE>"
    else:
        doc_str = ""
    return doc_str


def _safe_makedirs(path: Path) -> None:
    if not path.exists():
        os.makedirs(path)
        _logger.info("Created directory: %s", path)


_string_type_def: Final = TypeDef(
    instance_type="String",
    init="new PrimitiveLoader<String>(String.class)",
    name="StringInstance",
    loader_type="Loader<String>",
)

_int_type_def: Final = TypeDef(
    instance_type="Integer",
    init="new PrimitiveLoader<Integer>(Integer.class)",
    name="IntegerInstance",
    loader_type="Loader<Integer>",
)

_long_type_def: Final = TypeDef(
    instance_type="Long",
    name="LongInstance",
    loader_type="Loader<Long>",
    init="new PrimitiveLoader<Long>(Long.class)",
)

_float_type_def: Final = TypeDef(
    instance_type="Double",
    name="DoubleInstance",
    loader_type="Loader<Double>",
    init="new PrimitiveLoader<Double>(Double.class)",
)

_bool_type_def: Final = TypeDef(
    instance_type="Boolean",
    name="BooleanInstance",
    loader_type="Loader<Boolean>",
    init="new PrimitiveLoader<Boolean>(Boolean.class)",
)

_null_type_def: Final = TypeDef(
    instance_type="Object",
    name="NullInstance",
    loader_type="Loader<Object>",
    init="new NullLoader()",
)

_any_type_def: Final = TypeDef(
    instance_type="Object",
    name="AnyInstance",
    init="new AnyLoader()",
    loader_type="Loader<Object>",
)

prims: Final = {
    "http://www.w3.org/2001/XMLSchema#string": _string_type_def,
    "http://www.w3.org/2001/XMLSchema#int": _int_type_def,
    "http://www.w3.org/2001/XMLSchema#long": _long_type_def,
    "http://www.w3.org/2001/XMLSchema#float": _float_type_def,
    "http://www.w3.org/2001/XMLSchema#double": _float_type_def,
    "http://www.w3.org/2001/XMLSchema#boolean": _bool_type_def,
    "https://w3id.org/cwl/salad#null": _null_type_def,
    "https://w3id.org/cwl/salad#Any": _any_type_def,
    "string": _string_type_def,
    "int": _int_type_def,
    "long": _long_type_def,
    "float": _float_type_def,
    "double": _float_type_def,
    "boolean": _bool_type_def,
    "null": _null_type_def,
    "Any": _any_type_def,
}


class JavaCodeGen(CodeGenBase):
    def __init__(
        self,
        base: str,
        target: Optional[str],
        examples: Optional[str],
        package: str,
        copyright: Optional[str],
    ) -> None:
        super().__init__()
        self.base_uri: Final = base
        self.examples: Final = examples
        self.package: Final = package
        self.artifact: Final = self.package.split(".")[-1]
        self.copyright: Final = copyright
        self.target_dir: Final = Path(target or ".").resolve()
        rel_package_dir = self.package.replace(".", "/")
        self.rel_package_dir: Final = Path(rel_package_dir)
        self.main_src_dir: Final = self.target_dir / "src" / "main" / "java" / rel_package_dir
        self.test_src_dir: Final = self.target_dir / "src" / "test" / "java" / rel_package_dir
        self.test_resources_dir: Final = (
            self.target_dir / "src" / "test" / "resources" / rel_package_dir
        )

    def prologue(self) -> None:
        for src_dir in [self.main_src_dir, self.test_src_dir]:
            _safe_makedirs(src_dir)

        for primitive in prims.values():
            self.declare_type(primitive)

    @staticmethod
    def property_name(name: str) -> str:
        avn: Final = schema.avro_field_name(name)
        return avn

    @staticmethod
    def safe_name(name: str) -> str:
        avn = JavaCodeGen.property_name(name)
        if avn in ("class", "extends", "abstract", "default", "package"):
            # reserved words
            avn = avn + "_"
        if avn and avn.startswith("anon."):
            avn = avn[5:]
        return avn

    def interface_name(self, n: str) -> str:
        return self.safe_name(n)

    def begin_class(
        self,
        classname: str,
        extends: MutableSequence[str],
        doc: str,
        abstract: bool,
        field_names: MutableSequence[str],
        idfield: str,
        optional_fields: set[str],
    ) -> None:
        cls = self.interface_name(classname)
        self.current_class = cls
        self.current_class_is_abstract = abstract
        self.current_loader = StringIO()
        self.current_fieldtypes: dict[str, TypeDef] = {}
        self.current_fields = StringIO()
        interface_doc_str = f"* Auto-generated interface for <I>{classname}</I><BR>"
        if not abstract:
            implemented_by = "This interface is implemented by {{@link {}Impl}}<BR>"
            interface_doc_str += implemented_by.format(cls)
        interface_doc_str += _doc_to_doc_string(doc)
        class_doc_str = f"* Auto-generated class implementation for <I>{classname}</I><BR>"
        class_doc_str += _doc_to_doc_string(doc)
        target = self.main_src_dir / f"{cls}.java"
        with open(target, "w") as f:
            _logger.info("Writing file: %s", target)
            if extends:
                ext = "extends " + ", ".join(self.interface_name(e) for e in extends) + ", Saveable"
            else:
                ext = "extends Saveable"
            f.write(
                """// Copyright Common Workflow Language project contributors
"""
            )
            if self.copyright:
                f.write(
                    """// {copyright}
""".format(
                        copyright=self.copyright
                    )
                )
            f.write(
                """//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package {package};

import {package}.utils.LoadingOptions;
import {package}.utils.Saveable;

/**
{interface_doc_str}
 */
public interface {cls} {ext} {{

  java.util.Map<String, Object> getExtensionFields();
  LoadingOptions getLoadingOptions();
""".format(
                    package=self.package,
                    cls=cls,
                    ext=ext,
                    interface_doc_str=interface_doc_str,
                )
            )

        if self.current_class_is_abstract:
            return

        target = self.main_src_dir / f"{cls}Impl.java"
        with open(target, "w") as f:
            _logger.info("Writing file: %s", target)
            f.write(
                """// Copyright Common Workflow Language project contributors
"""
            )
            if self.copyright:
                f.write(
                    """// {copyright}
""".format(
                        copyright=self.copyright
                    )
                )
            f.write(
                """//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package {package};

import {package}.utils.LoaderInstances;
import {package}.utils.LoadingOptions;
import {package}.utils.LoadingOptionsBuilder;
import {package}.utils.SaveableImpl;
import {package}.utils.ValidationException;

/**
{class_doc_str}
 */
public class {cls}Impl extends SaveableImpl implements {cls} {{
  private LoadingOptions loadingOptions_ = new LoadingOptionsBuilder().build();
  private java.util.Map<String, Object> extensionFields_ =
      new java.util.HashMap<String, Object>();
  public LoadingOptions getLoadingOptions() {{
    return this.loadingOptions_;
  }}
  public java.util.Map<String, Object> getExtensionFields() {{
    return this.extensionFields_;
  }}
""".format(
                    package=self.package,
                    cls=cls,
                    class_doc_str=class_doc_str,
                )
            )
        self.current_loader.write(
            """
  /**
   * Used by {{@link {package}.utils.RootLoader}} to construct instances of {cls}Impl.
   *
   * @param __doc_            Document fragment to load this record object from (presumably a
                              {{@link java.util.Map}}).
   * @param __baseUri_        Base URI to generate child document IDs against.
   * @param __loadingOptions  Context for loading URIs and populating objects.
   * @param __docRoot_        ID at this position in the document (if available) (maybe?)
   * @throws ValidationException If the document fragment is not a {{@link java.util.Map}}
   *                             or validation of fields fails.
   */
  public {cls}Impl(
      final Object __doc_,
      final String __baseUri_,
      LoadingOptions __loadingOptions,
      final String __docRoot_) {{
    super(__doc_, __baseUri_, __loadingOptions, __docRoot_);
    // Prefix plumbing variables with '__' to reduce likelihood of collision with
    // generated names.
    String __baseUri = __baseUri_;
    String __docRoot = __docRoot_;
    if (!(__doc_ instanceof java.util.Map)) {{
      throw new ValidationException("{cls}Impl called on non-map");
    }}
    final java.util.Map<String, Object> __doc = (java.util.Map<String, Object>) __doc_;
    final java.util.List<ValidationException> __errors =
        new java.util.ArrayList<ValidationException>();
    if (__loadingOptions != null) {{
      this.loadingOptions_ = __loadingOptions;
    }}
""".format(
                cls=cls, package=self.package
            )
        )

    def end_class(self, classname: str, field_names: list[str]) -> None:
        """Finish this class."""
        with open(self.main_src_dir / f"{self.current_class}.java", "a") as f:
            f.write(
                """
}
"""
            )
        if self.current_class_is_abstract:
            return

        self.current_loader.write(
            """    if (!__errors.isEmpty()) {
      throw new ValidationException("Trying 'RecordField'", __errors);
    }
"""
        )
        for fieldname in field_names:
            fieldtype = self.current_fieldtypes.get(fieldname)
            if fieldtype is None:
                continue
            self.current_loader.write(
                """    this.{safename} = ({type}) {safename};
""".format(
                    safename=self.safe_name(fieldname), type=fieldtype.instance_type
                )
            )
        self.current_loader.write(
            """    for (String field:__doc.keySet()) {
      if (!attrs.contains(field)) {
        if (field.contains(":")) {
          String expanded_field = __loadingOptions.expandUrl(field, "", false, false, null);
          extensionFields_.put(expanded_field, __doc.get(field));
        }
      }
    }
"""
        )
        self.current_loader.write("""  }""")
        target = self.main_src_dir / f"{self.current_class}Impl.java"
        with open(
            target,
            "a",
        ) as f:
            f.write(self.current_fields.getvalue())
            f.write(self.current_loader.getvalue())
            f.write(
                f"""
  private java.util.List<String> attrs = java.util.Arrays.asList("{'", "'.join(field_names)}");
}}
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
            sub = [self.type_loader(i) for i in type_declaration]
            if len(sub) < 2:
                return sub[0]

            if len(sub) == 2:
                type_1 = sub[0]
                type_2 = sub[1]
                type_1_name = type_1.name
                type_2_name = type_2.name
                if type_1_name == "NullInstance" or type_2_name == "NullInstance":
                    non_null_type = type_1 if type_1.name != "NullInstance" else type_2
                    return self.declare_type(
                        TypeDef(
                            instance_type="java.util.Optional<{}>".format(
                                non_null_type.instance_type
                            ),
                            init=f"new OptionalLoader({non_null_type.name})",
                            name=f"optional_{non_null_type.name}",
                            loader_type="Loader<java.util.Optional<{}>>".format(
                                non_null_type.instance_type
                            ),
                        )
                    )
                if (
                    type_1_name == f"array_of_{type_2_name}"
                    or type_2_name == f"array_of_{type_1_name}"
                ) and USE_ONE_OR_LIST_OF_TYPES:
                    if type_1_name == f"array_of_{type_2_name}":
                        single_type = type_2
                        array_type = type_1
                    else:
                        single_type = type_1
                        array_type = type_2
                    fqclass = f"{self.package}.{single_type.instance_type}"
                    return self.declare_type(
                        TypeDef(
                            instance_type=f"{self.package}.utils.OneOrListOf<{fqclass}>",
                            init="new OneOrListOfLoader<{}>({}, {})".format(
                                fqclass, single_type.name, array_type.name
                            ),
                            name=f"one_or_array_of_{single_type.name}",
                            loader_type="Loader<{}.utils.OneOrListOf<{}>>".format(
                                self.package, fqclass
                            ),
                        )
                    )
            return self.declare_type(
                TypeDef(
                    instance_type="Object",
                    init="new UnionLoader(new Loader[] {{ {} }})".format(
                        ", ".join(s.name for s in sub)
                    ),
                    name="union_of_{}".format("_or_".join(s.name for s in sub)),
                    loader_type="Loader<Object>",
                )
            )
        if isinstance(type_declaration, MutableMapping):
            if type_declaration["type"] in (
                "array",
                "https://w3id.org/cwl/salad#array",
            ):
                i = self.type_loader(type_declaration["items"])
                instance_type = (
                    "java.util.List<String>"
                    if i.instance_type == "String"
                    else "java.util.List<Object>"
                )
                return self.declare_type(
                    TypeDef(
                        # special doesn't work out with subclassing, gotta be more clever
                        # instance_type="List<{}>".format(i.instance_type),
                        instance_type=instance_type,
                        name=f"array_of_{i.name}",
                        init=f"new ArrayLoader({i.name})",
                        loader_type=f"Loader<java.util.List<{i.instance_type}>>",
                    )
                )
            if type_declaration["type"] in (
                "map",
                "https://w3id.org/cwl/salad#map",
            ):
                i = self.type_loader(type_declaration["values"])
                return self.declare_type(
                    TypeDef(
                        # special doesn't work out with subclassing, gotta be more clever
                        # instance_type="Map<String, {}>".format(i.instance_type),
                        instance_type=f"java.util.Map<String, {i.instance_type}>",
                        name=f"map_of_{i.name}",
                        init="new MapLoader({}, {}, {})".format(
                            i.name,
                            (
                                f'"{container}"' if container is not None else self.to_java(None)
                            ),  # noqa: B907
                            self.to_java(no_link_check),
                        ),
                        loader_type=f"Loader<java.util.Map<String, {i.instance_type}>>",
                    )
                )
            if type_declaration["type"] in ("enum", "https://w3id.org/cwl/salad#enum"):
                return self.type_loader_enum(type_declaration)
            if type_declaration["type"] in (
                "record",
                "https://w3id.org/cwl/salad#record",
            ):
                is_abstract = type_declaration.get("abstract", False)
                fqclass = "{}.{}".format(self.package, self.safe_name(type_declaration["name"]))
                return self.declare_type(
                    TypeDef(
                        instance_type=self.safe_name(type_declaration["name"]),
                        name=self.safe_name(type_declaration["name"]),
                        init="new RecordLoader<{clazz}>({clazz}{ext}.class, "
                        "{container}, {no_link_check})".format(
                            clazz=fqclass,
                            ext="Impl" if not is_abstract else "",
                            container=(
                                f'"{container}"'
                                if container is not None
                                else self.to_java(None)  # noqa: B907
                            ),
                            no_link_check=self.to_java(no_link_check),
                        ),
                        loader_type=f"Loader<{fqclass}>",
                    )
                )
            if type_declaration["type"] in (
                "union",
                "https://w3id.org/cwl/salad#union",
            ):
                # Declare the named loader to handle recursive union definitions
                loader_name = self.safe_name(type_declaration["name"])
                loader_type = TypeDef(
                    instance_type="Object",
                    init="new UnionLoader(new Loader[] {})",
                    name=loader_name,
                    loader_type="Loader<Object>",
                )
                self.declare_type(loader_type)
                # Parse inner types
                sub = [self.type_loader(i) for i in type_declaration["names"]]

                if len(sub) == 2:
                    type_1 = sub[0]
                    type_2 = sub[1]
                    type_1_name = type_1.name
                    type_2_name = type_2.name
                    if type_1_name == "NullInstance" or type_2_name == "NullInstance":
                        non_null_type = type_1 if type_1.name != "NullInstance" else type_2
                        sub = [
                            TypeDef(
                                instance_type="java.util.Optional<{}>".format(
                                    non_null_type.instance_type
                                ),
                                init=f"new OptionalLoader({non_null_type.name})",
                                name=f"optional_{non_null_type.name}",
                                loader_type="Loader<java.util.Optional<{}>>".format(
                                    non_null_type.instance_type
                                ),
                            )
                        ]
                    elif (
                        type_1_name == f"array_of_{type_2_name}"
                        or type_2_name == f"array_of_{type_1_name}"
                    ) and USE_ONE_OR_LIST_OF_TYPES:
                        if type_1_name == f"array_of_{type_2_name}":
                            single_type = type_2
                            array_type = type_1
                        else:
                            single_type = type_1
                            array_type = type_2
                        fqclass = f"{self.package}.{single_type.instance_type}"
                        sub = [
                            TypeDef(
                                instance_type=f"{self.package}.utils.OneOrListOf<{fqclass}>",
                                init="new OneOrListOfLoader<{}>({}, {})".format(
                                    fqclass, single_type.name, array_type.name
                                ),
                                name=f"one_or_array_of_{single_type.name}",
                                loader_type="Loader<{}.utils.OneOrListOf<{}>>".format(
                                    self.package, fqclass
                                ),
                            )
                        ]
                # Register lazy initialization for the loader
                self.add_lazy_init(
                    LazyInitDef(
                        loader_name,
                        "((UnionLoader) {}).addLoaders(new Loader[] {{ {} }});".format(
                            loader_name, ", ".join(s.name for s in sub)
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
                    name=self.safe_name(type_declaration) + "Loader",
                    init="new ExpressionLoader()",
                    loader_type="Loader<String>",
                    instance_type="String",
                )
            )
        return self.collected_types[self.safe_name(type_declaration)]

    def type_loader_enum(self, type_declaration: dict[str, Any]) -> TypeDef:
        """Build an enum type loader for the given declaration."""
        symbols = [self.property_name(sym) for sym in type_declaration["symbols"]]
        for sym in symbols:
            self.add_vocab(shortname(sym), sym)
        clazz = self.safe_name(type_declaration["name"])
        symbols_decl = 'new String[] {{"{}"}}'.format('", "'.join(sym for sym in symbols))
        enum_path = self.main_src_dir / f"{clazz}.java"
        with open(enum_path, "w") as f:
            _logger.info("Writing file: %s", enum_path)
            f.write(
                """// Copyright Common Workflow Language project contributors
"""
            )
            if self.copyright:
                f.write(
                    """// {copyright}
""".format(
                        copyright=self.copyright
                    )
                )
            f.write(
                """//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//      http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package {package};

import {package}.utils.ValidationException;

public enum {clazz} {{
""".format(
                    package=self.package, clazz=clazz
                )
            )
            for i, sym in enumerate(symbols):
                suffix = "," if i < (len(symbols) - 1) else ";"
                const = self.safe_name(sym).replace("-", "_").replace(".", "_").upper()
                f.write(f"""  {const}("{sym}"){suffix}\n""")  # noqa: B907
            f.write(
                """
  private static String[] symbols = {symbols_decl};
  private String docVal;

  private {clazz}(final String docVal) {{
    this.docVal = docVal;
  }}

  public static {clazz} fromDocumentVal(final String docVal) {{
    for(final {clazz} val : {clazz}.values()) {{
      if(val.docVal.equals(docVal)) {{
        return val;
      }}
    }}
    throw new ValidationException(String.format("Expected one of %s", {clazz}.symbols, docVal));
  }}
}}
""".format(
                    clazz=clazz, symbols_decl=symbols_decl
                )
            )
        return self.declare_type(
            TypeDef(
                instance_type=clazz,
                name=self.safe_name(type_declaration["name"]),
                loader_type=f"Loader<{clazz}>",
                init=f"new EnumLoader({clazz}.class)",
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
        fieldname = name
        property_name = self.property_name(fieldname)
        cap_case_property_name = property_name[0].upper() + property_name[1:]
        if cap_case_property_name == "Class":
            cap_case_property_name = "Class_"

        safename = self.safe_name(fieldname)
        self.current_fieldtypes[property_name] = fieldtype
        getter_doc_str = """  /**
   * Getter for property <I>{fieldname}</I><BR>
{field_doc_str}
   */
""".format(
            fieldname=fieldname, field_doc_str=_doc_to_doc_string(doc, indent_level=1)
        )
        target = self.main_src_dir / f"{self.current_class}.java"
        with open(target, "a") as f:
            f.write(
                """
{getter_doc_str}
  {type} get{capfieldname}();""".format(
                    getter_doc_str=getter_doc_str,
                    capfieldname=cap_case_property_name,
                    type=fieldtype.instance_type,
                )
            )

        if self.current_class_is_abstract:
            return

        self.current_fields.write(
            """
  private {type} {safename};

{getter_doc_str}
  public {type} get{capfieldname}() {{
    return this.{safename};
  }}
""".format(
                safename=safename,
                capfieldname=cap_case_property_name,
                getter_doc_str=getter_doc_str,
                type=fieldtype.instance_type,
            )
        )

        self.current_loader.write(
            """    {type} {safename};
""".format(
                type=fieldtype.instance_type, safename=safename
            )
        )
        if optional:
            self.current_loader.write(
                """
    if (__doc.containsKey("{fieldname}")) {{
""".format(
                    fieldname=property_name
                )
            )
            spc = "  "
        else:
            spc = ""

        self.current_loader.write(
            """{spc}    try {{
{spc}      {safename} =
{spc}          LoaderInstances
{spc}              .{fieldtype}
{spc}              .loadField(__doc.get("{fieldname}"), __baseUri, __loadingOptions);
{spc}    }} catch (ValidationException e) {{
{spc}      {safename} = null; // won't be used but prevents compiler from complaining.
{spc}      final String __message = "the `{fieldname}` field is not valid because:";
{spc}      __errors.add(new ValidationException(__message, e));
{spc}    }}
""".format(
                fieldtype=fieldtype.name,
                safename=safename,
                fieldname=property_name,
                spc=spc,
            )
        )

        if optional:
            self.current_loader.write(
                """
    }} else {{
      {safename} = null;
    }}
""".format(
                    safename=safename
                )
            )

    def declare_id_field(
        self,
        name: str,
        fieldtype: TypeDef,
        doc: Optional[str],
        optional: bool,
    ) -> None:
        if self.current_class_is_abstract:
            return

        self.declare_field(name, fieldtype, doc, True, "")
        if optional:
            set_uri = """
    Boolean __original_is_null = {safename} == null;
    if ({safename} == null) {{
      if (__docRoot != null) {{
        {safename} = java.util.Optional.of(__docRoot);
      }} else {{
        {safename} = java.util.Optional.of("_:" + java.util.UUID.randomUUID().toString());
      }}
    }}
    if (__original_is_null) {{
        __baseUri = __baseUri_;
    }} else {{
        __baseUri = (String) {safename}.orElse(null);
    }}
"""
        else:
            set_uri = """
    if ({safename} == null) {{
      if (__docRoot != null) {{
        {safename} = __docRoot;
      }} else {{
        throw new ValidationException("Missing {fieldname}");
      }}
    }}
    __baseUri = (String) {safename};
"""

        self.current_loader.write(
            set_uri.format(safename=self.safe_name(name), fieldname=shortname(name))
        )

    def uri_loader(
        self,
        inner: TypeDef,
        scoped_id: bool,
        vocab_term: bool,
        ref_scope: Optional[int],
        no_link_check: Optional[bool] = None,
    ) -> TypeDef:
        instance_type = inner.instance_type or "Object"
        return self.declare_type(
            TypeDef(
                instance_type=instance_type,  # ?
                name=f"uri_{inner.name}_{scoped_id}_{vocab_term}_{ref_scope}_{no_link_check}",
                init="new UriLoader({}, {}, {}, {}, {})".format(
                    inner.name,
                    self.to_java(scoped_id),
                    self.to_java(vocab_term),
                    self.to_java(ref_scope),
                    self.to_java(no_link_check),
                ),
                is_uri=True,
                scoped_id=scoped_id,
                ref_scope=ref_scope,
                loader_type=f"Loader<{instance_type}>",
            )
        )

    def idmap_loader(
        self, field: str, inner: TypeDef, map_subject: str, map_predicate: Optional[str]
    ) -> TypeDef:
        instance_type: Final = inner.instance_type or "Object"
        return self.declare_type(
            TypeDef(
                instance_type=instance_type,
                name=f"idmap_{self.safe_name(field)}_{inner.name}",
                init='new IdMapLoader({}, "{}", "{}")'.format(
                    inner.name, map_subject, map_predicate
                ),
                loader_type=f"Loader<{instance_type}>",
            )
        )

    def typedsl_loader(self, inner: TypeDef, ref_scope: Union[int, None]) -> TypeDef:
        """Construct the TypeDef for the given DSL loader."""
        instance_type = inner.instance_type or "Object"
        return self.declare_type(
            TypeDef(
                instance_type=instance_type,
                name=f"typedsl_{inner.name}_{ref_scope}",
                init=f"new TypeDslLoader({inner.name}, {ref_scope})",
                loader_type=f"Loader<{instance_type}>",
            )
        )

    def to_java(self, val: Any) -> Any:
        if val is True:
            return "true"
        elif val is None:
            return "null"
        elif val is False:
            return "false"
        return val

    def epilogue(self, root_loader: TypeDef) -> None:
        pd = "This project contains Java objects and utilities "
        pd = pd + ' auto-generated by <a href="https://github.com/'
        pd = pd + 'common-workflow-language/schema_salad">Schema Salad</a>'
        pd = pd + " for parsing documents corresponding to the "
        pd = pd + str(self.base_uri) + " schema."

        template_vars: MutableMapping[str, str] = dict(
            base_uri=self.base_uri,
            package=self.package,
            group_id=self.package,
            artifact_id=self.artifact,
            version="0.0.1-SNAPSHOT",
            project_name=self.package,
            project_description=pd,
            license_name="Apache License, Version 2.0",
            license_url="https://www.apache.org/licenses/LICENSE-2.0.txt",
        )

        def template_from_resource(resource: Traversable) -> string.Template:
            template_str: Final = resource.read_text("utf-8")
            template: Final = string.Template(template_str)
            return template

        def expand_resource_template_to(resource: str, path: Path) -> None:
            template: Final = template_from_resource(
                files("schema_salad").joinpath(f"java/{resource}")
            )
            src: Final = template.safe_substitute(template_vars)
            _ensure_directory_and_write(path, src)

        expand_resource_template_to("pom.xml", self.target_dir / "pom.xml")
        expand_resource_template_to("gitignore", self.target_dir / ".gitignore")
        expand_resource_template_to("package.html", self.main_src_dir / "package.html")
        expand_resource_template_to(
            "overview.html",
            self.target_dir / "src" / "main" / "javadoc" / "overview.html",
        )
        expand_resource_template_to(
            "MANIFEST.MF",
            self.target_dir / "src" / "main" / "resources" / "META-INF" / "MANIFEST.MF",
        )
        expand_resource_template_to("README.md", self.target_dir / "README.md")

        vocab = ""
        rvocab = ""
        for k in sorted(self.vocab.keys()):
            vocab += f"""    vocab.put("{k}", "{self.vocab[k]}");\n"""  # noqa: B907
            rvocab += f"""    rvocab.put("{self.vocab[k]}", "{k}");\n"""  # noqa: B907

        loader_instances = ""
        for _, collected_type in self.collected_types.items():
            loader_instances += "  public static {} {} = {};\n".format(
                collected_type.loader_type, collected_type.name, collected_type.init
            )

        if self.lazy_inits:
            loader_instances += "\n  static {\n"
            for lazy_init in self.lazy_inits.values():
                loader_instances += f"    {lazy_init.init}\n"
            loader_instances += "  }\n"

        example_tests = ""
        if self.examples:
            _safe_makedirs(self.test_resources_dir)
            utils_resources = self.test_resources_dir / "utils"
            if os.path.exists(utils_resources):
                shutil.rmtree(utils_resources)
            shutil.copytree(self.examples, utils_resources)
            for example_name in os.listdir(self.examples):
                if example_name.startswith("valid"):
                    basename = os.path.basename(example_name).rsplit(".", 1)[0]
                    basename = re.sub(BASIC_JAVA_IDENTIFIER_RE, "_", basename)
                    example_tests += """
  @org.junit.Test
  public void test{basename}ByString() throws Exception {{
    java.net.URL url = getClass().getResource("{example_name}");
    java.nio.file.Path resPath = java.nio.file.Paths.get(url.toURI());
    String yaml = new String(java.nio.file.Files.readAllBytes(resPath), "UTF8");
    RootLoader.loadDocument(yaml, url.toString());
  }}

  @org.junit.Test
  public void test{basename}ByPath() throws Exception {{
    java.net.URL url = getClass().getResource("{example_name}");
    java.nio.file.Path resPath = java.nio.file.Paths.get(url.toURI());
    RootLoader.loadDocument(resPath);
  }}

  @org.junit.Test
  public void test{basename}ByMap() throws Exception {{
    java.net.URL url = getClass().getResource("{example_name}");
    java.nio.file.Path resPath = java.nio.file.Paths.get(url.toURI());
    String yaml = new String(java.nio.file.Files.readAllBytes(resPath), "UTF8");
    java.util.Map<String, Object> doc;
    doc = (java.util.Map<String, Object>) YamlUtils.mapFromString(yaml);
    RootLoader.loadDocument(doc, url.toString());
  }}""".format(
                        basename=basename,
                        example_name=example_name,
                    )

        template_args: MutableMapping[str, str] = dict(
            package=self.package,
            vocab=vocab,
            rvocab=rvocab,
            loader_instances=loader_instances,
            root_loader_name=root_loader.name,
            root_loader_instance_type=root_loader.instance_type or "Object",
            example_tests=example_tests,
        )

        util_src_dirs: Final = {
            "main_utils": self.main_src_dir,
            "test_utils": self.test_src_dir,
        }
        for util_src, util_target in util_src_dirs.items():
            for util in files("schema_salad").joinpath(f"java/{util_src}").iterdir():
                src_path = util_target / "utils" / util.name
                src_template = template_from_resource(util)
                src = src_template.safe_substitute(template_args)
                _ensure_directory_and_write(src_path, src)

    def secondaryfilesdsl_loader(self, inner: TypeDef) -> TypeDef:
        instance_type: Final = inner.instance_type or "Object"
        return self.declare_type(
            TypeDef(
                instance_type=instance_type,
                name=f"secondaryfilesdsl_{inner.name}",
                init=f"new SecondaryFilesDslLoader({inner.name})",
                loader_type=f"Loader<{instance_type}>",
            )
        )
