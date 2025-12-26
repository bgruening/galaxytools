"""
C++17 code generator for a given Schema Salad definition.

Currently only supports emitting YAML from the C++ objects, not yet parsing
YAML into C++ objects.

The generated code requires the libyaml-cpp library & headers

To see an example of usage, look at schema_salad/tests/codegen/cwl.cpp
which can be combined with the CWL V1.0 schema as shown below::

  schema-salad-tool --codegen cpp \
          schema_salad/tests/test_schema/CommonWorkflowLanguage.yml \
          > cwl_v1_0.h

  g++ --std=c++20 -I. -lyaml-cpp schema_salad/tests/codegen/cwl.cpp -o cwl-v1_0-test
  ./cwl-v1_0-test

  # g++ versions older than version 10 may need "--std=c++2a" instead of "--std=c++20"
"""

import os
import re
from typing import IO, Any, Optional, Union, cast

from . import _logger
from .codegen_base import CodeGenBase, TypeDef
from .exceptions import SchemaException
from .schema import shortname
from .utils import aslist


def q(s: str) -> str:
    """Put quotes around a string."""
    return '"' + s + '"'


def replaceKeywords(s: str) -> str:
    """Rename keywords that are reserved in C++."""
    if s in (
        "class",
        "enum",
        "int",
        "long",
        "float",
        "double",
        "default",
        "stdin",
        "stdout",
        "stderr",
        "union",
    ):
        s = s + "_"
    return s


def safename(name: str) -> str:
    """Create a C++ safe name."""
    classname = re.sub(r"[^a-zA-Z0-9]", "_", name)
    return replaceKeywords(classname)


def safenamespacename(name: str) -> str:
    """Create a C++ safe name for namespaces."""
    name = re.sub(r"^[a-zA-Z0-9]+://", "", name)  # remove protocol
    name = re.sub(r"//+", "", name)  # remove duplicate slashes
    name = re.sub(r"/$", "", name)  # remove trailing slashes
    name = re.sub(r"[^a-zA-Z0-9/]", "_", name)
    name = re.sub(r"[/]", "::", name)

    return name


# TODO: this should be somehow not really exists
def safename2(name: dict[str, str]) -> str:
    """Create a namespaced safename."""
    return safenamespacename(name["namespace"]) + "::" + safename(name["classname"])


def split_name(s: str) -> tuple[str, str]:
    """Split url name into its components.

    Splits names like https://xyz.xyz/blub#cwl/class
    into its class path and non class path
    """
    t = s.split("#")
    if len(t) != 2:
        raise ValueError("Expected field to be formatted as 'https://xyz.xyz/blub#cwl/class'.")
    return (t[0], t[1])


def split_field(s: str) -> tuple[str, str, str]:
    """Split field into its components.

    similar to split_name but for field names
    """
    (namespace, field) = split_name(s)
    t = field.split("/")
    if len(t) != 2:
        raise ValueError("Expected field to be formatted as 'https://xyz.xyz/blub#cwl/class'.")
    return (namespace, t[0], t[1])


class ClassDefinition:
    """Prototype of a class."""

    def __init__(self, name: str):
        """Initialize the class definition with a name."""
        self.fullName = name
        self.extends: list[dict[str, str]] = []

        # List of types from parent classes that have been specialized
        self.specializationTypes: list[str] = []

        # this includes fields that are also inheritant
        self.allfields: list[FieldDefinition] = []
        self.fields: list[FieldDefinition] = []
        self.abstract = False
        (self.namespace, self.classname) = split_name(name)
        self.namespace = safenamespacename(self.namespace)
        self.classname = safename(self.classname)

    def writeFwdDeclaration(self, target: IO[str], fullInd: str, ind: str) -> None:
        """Write forward declaration."""
        target.write(f"{fullInd}namespace {self.namespace} {{ struct {self.classname}; }}\n")

    def writeDefinition(
        self, target: IO[Any], fullInd: str, ind: str, common_namespace: str
    ) -> None:
        """Write definition of the class."""
        target.write(f"{fullInd}namespace {self.namespace} {{\n")
        target.write(f"{fullInd}struct {self.classname}")
        extends = list(map(safename2, self.extends))
        override = ""
        virtual = "virtual "
        if len(self.extends) > 0:
            target.write(f"\n{fullInd}{ind}: ")
            target.write(f"\n{fullInd}{ind}, ".join(extends))
            override = " override"
            virtual = ""
        target.write(" {\n")

        for field in self.fields:
            field.writeDefinition(target, fullInd + ind, ind, self.namespace)

        if self.abstract:
            target.write(f"{fullInd}{ind}virtual ~{self.classname}() = 0;\n")
        else:
            target.write(f"{fullInd}{ind}{virtual}~{self.classname}(){override} = default;\n")

        target.write(
            f"{fullInd}{ind}{virtual}auto toYaml([[maybe_unused]] "
            f"{common_namespace}::store_config const& config) const -> YAML::Node{override};\n"
        )
        target.write(f"{fullInd}{ind}{virtual}void fromYaml(YAML::Node const& n){override};\n")
        target.write(f"{fullInd}}};\n")
        target.write(f"{fullInd}}}\n\n")

    def writeImplDefinition(
        self, target: IO[str], fullInd: str, ind: str, common_namespace: str
    ) -> None:
        """Write definition with implementation."""
        extends = list(map(safename2, self.extends))

        # Declaring default destructor
        if self.abstract:
            target.write(
                f"{fullInd}inline {self.namespace}::{self.classname}::~{self.classname}() = default;\n"
            )

        # Write toYaml function
        target.write(
            f"{fullInd}inline auto {self.namespace}::{self.classname}::toYaml([[maybe_unused]] "
            f"::{common_namespace}::store_config const& config) const -> YAML::Node {{\n"
            f"{fullInd}{ind}using ::{common_namespace}::toYaml;\n"
            f"{fullInd}{ind}auto n = YAML::Node{{}};\n"
            f"{fullInd}{ind}if (config.generateTags) {{\n"
            f'{fullInd}{ind}{ind}n.SetTag("{self.classname}");\n'
            f"{fullInd}{ind}}}\n"
        )
        for e in extends:
            target.write(f"{fullInd}{ind}n = mergeYaml(n, {e}::toYaml(config));\n")

        for field in self.fields:
            fieldname = safename(field.name)

            target.write(f"{fullInd}{ind}{{\n")
            target.write(f"{fullInd}{ind}{ind} auto member = toYaml(*{fieldname}, config);\n")
            if field.typeDSL:
                target.write(f"{fullInd}{ind}{ind} member = simplifyType(member, config);\n")
            target.write(
                f"{fullInd}{ind}{ind} member = convertListToMap(member, "
                f"{q(field.mapSubject)}, {q(field.mapPredicate)}, config);\n"
            )
            target.write(f"{fullInd}{ind}{ind}addYamlField(n, {q(field.name)}, member);\n")
            target.write(f"{fullInd}{ind}}}\n")

        target.write(f"{fullInd}{ind}return n;\n{fullInd}}}\n")

        # Write fromYaml function
        functionname = f"{self.namespace}::{self.classname}::fromYaml"
        target.write(
            f"{fullInd}inline void {functionname}([[maybe_unused]] YAML::Node const& n) {{\n"
            f"{fullInd}{ind}using ::{common_namespace}::fromYaml;\n"
        )
        for e in extends:
            target.write(f"{fullInd}{ind}{e}::fromYaml(n);\n")

        for field in self.fields:
            fieldname = safename(field.name)
            expandType = ""
            if field.typeDSL:
                expandType = "expandType"

            target.write(
                f"{fullInd}{ind}{{\n"
                f"{fullInd}{ind}{ind}auto nodeAsList = convertMapToList("
                f"n[{q(field.name)}], {q(field.mapSubject)}, {q(field.mapPredicate)});\n"
                f"{fullInd}{ind}{ind}auto expandedNode = {expandType}(nodeAsList);\n"
                f"{fullInd}{ind}{ind}fromYaml(expandedNode, *{fieldname});\n"
                f"{fullInd}{ind}}}\n"
            )

        target.write(f"{fullInd}}}\n")

        # write type detection function
        if not self.abstract:
            e = f"::{self.namespace}::{self.classname}"
            target.write(
                f"namespace {common_namespace} {{\n"
                f"template <>\n"
                f"struct DetectAndExtractFromYaml<{e}> {{\n"
                f"    auto operator()(YAML::Node const& n) const -> std::optional<{e}> {{\n"
                f"        if (!n.IsDefined()) return std::nullopt;\n"
                f"        if (!n.IsMap()) return std::nullopt;\n"
                f"        auto res = {e}{{}};\n\n"
            )
            for field in self.fields:
                fieldname = safename(field.name)
                target.write(
                    f"        if constexpr (::{common_namespace}::IsConstant<"
                    f"decltype(res.{fieldname})::value_t>::value) try {{\n"
                    f"            fromYaml(n[{q(field.name)}], *res.{fieldname});\n"
                    "            fromYaml(n, res);\n"
                    "            return res;\n"
                    f"        }} catch(...) {{}}\n\n"
                )
            target.write("        return std::nullopt;\n    }\n};\n}\n")


class FieldDefinition:
    """Prototype of a single field from a class definition."""

    def __init__(
        self,
        name: str,
        typeStr: str,
        optional: bool,
        mapSubject: str,
        mapPredicate: str,
        typeDSL: bool,
    ):
        """Initialize field definition.

        Creates a new field with name, its type, optional and which field to use to convert
        from list to map (or empty if it is not possible)
        """
        self.name = name
        self.typeStr = typeStr
        self.optional = optional
        self.mapSubject = mapSubject
        self.mapPredicate = mapPredicate
        self.typeDSL = typeDSL

    def writeDefinition(self, target: IO[Any], fullInd: str, ind: str, namespace: str) -> None:
        """Write a C++ definition for the class field."""
        name = safename(self.name)
        typeStr = self.typeStr.replace(namespace + "::", "")
        target.write(f"{fullInd}heap_object<{typeStr}> {name};\n")


class MapDefinition:
    """Prototype of a map."""

    def __init__(self, name: str, values: list[str]):
        """Initialize union definition with a name and possible values."""
        self.values = values
        (self.namespace, self.classname) = split_name(name)
        self.namespace = safenamespacename(self.namespace)
        self.classname = safename(self.classname)

    def _remove_namespace(self, typeStr: str) -> str:
        return typeStr.replace(f"{self.namespace}::", "")

    def writeFwdDeclaration(self, target: IO[str], fullInd: str, ind: str) -> None:
        """Write forward declaration."""
        target.write(f"{fullInd}namespace {self.namespace} {{ struct {self.classname}; }}\n")

    def writeDefinition(self, target: IO[str], ind: str, common_namespace: str) -> None:
        """Write map definition to output."""
        target.write(f"namespace {self.namespace} {{\n")
        if len(self.values) == 1:
            valueType = self._remove_namespace(self.values[0])
        else:
            valueType = f"std::variant<{', '.join(self._remove_namespace(v) for v in self.values)}>"
        target.write(f"struct {self.classname} {{\n")
        target.write(f"{ind}heap_object<std::map<std::string, {valueType}>> value;\n")
        target.write(
            f"{ind}auto toYaml([[maybe_unused]] "
            f"::{common_namespace}::store_config const& config) const -> YAML::Node;\n"
        )
        target.write(f"{ind}void fromYaml(YAML::Node const& n);\n")
        target.write("};\n")
        target.write("}\n\n")

    def writeImplDefinition(
        self, target: IO[str], fullInd: str, ind: str, common_namespace: str
    ) -> None:
        """Write definition with implementation."""
        # Write toYaml function
        functionname = f"{self.namespace}::{self.classname}::toYaml"
        target.write(
            f"{fullInd}inline auto {functionname}([[maybe_unused]] "
            f"::{common_namespace}::store_config const& config) const -> YAML::Node {{\n"
            f"{fullInd}{ind}using ::{common_namespace}::toYaml;\n"
            f"{fullInd}{ind}return toYaml(*value, config);\n"
            f"{fullInd}}}\n"
        )

        # Write fromYaml function
        functionname = f"{self.namespace}::{self.classname}::fromYaml"
        target.write(
            f"{fullInd}inline void {functionname}([[maybe_unused]] YAML::Node const& n) {{\n"
            f"{fullInd}{ind}using ::{common_namespace}::fromYaml;\n"
            f"{fullInd}{ind}fromYaml(n, *value);\n"
            f"{fullInd}}}\n"
        )


class UnionDefinition:
    """Prototype of a union."""

    def __init__(self, name: str, types: list[str]):
        """Initialize union definition with a name and possible types."""
        (self.namespace, self.classname) = split_name(name)
        self.namespace = safenamespacename(self.namespace)
        self.classname = safename(self.classname)
        self.types = (
            self._remove_namespace(types[0])
            if len(types) == 1
            else f"std::variant<{', '.join(self._remove_namespace(t) for t in types)}>"
        )

    def _remove_namespace(self, typeStr: str) -> str:
        return typeStr.replace(f"{self.namespace}::", "")

    def writeFwdDeclaration(self, target: IO[str], fullInd: str, ind: str) -> None:
        """Write forward declaration."""
        target.write(f"{fullInd}namespace {self.namespace} {{ struct {self.classname}; }}\n")

    def writeDefinition(self, target: IO[str], ind: str, common_namespace: str) -> None:
        """Write union definition to output."""
        target.write(f"namespace {self.namespace} {{\n")
        target.write(f"struct {self.classname} {{\n")
        target.write(f"{ind}{self.types} *value = nullptr;\n")
        target.write(f"{ind}{self.classname}();\n")
        target.write(f"{ind}~{self.classname}();\n")
        target.write(
            f"{ind}auto toYaml([[maybe_unused]] "
            f"::{common_namespace}::store_config const& config) const -> YAML::Node;\n"
        )
        target.write(f"{ind}void fromYaml(YAML::Node const& n);\n")
        target.write("};\n")
        target.write("}\n\n")

    def writeImplDefinition(
        self, target: IO[str], fullInd: str, ind: str, common_namespace: str
    ) -> None:
        """Write definition with implementation."""
        # Write constructor
        functionname = f"{self.namespace}::{self.classname}::{self.classname}"
        target.write(
            f"{fullInd}{functionname}() {{\n"
            f"{fullInd}{ind}value = new {self.types}();\n"
            f"{fullInd}}}\n"
        )

        # Write destructor
        functionname = f"{self.namespace}::{self.classname}::~{self.classname}"
        target.write(
            f"{fullInd}{functionname}() {{\n"
            f"{fullInd}{ind}if (value != nullptr) {{\n"
            f"{fullInd}{ind}{ind}delete value;\n"
            f"{fullInd}{ind}{ind}value = nullptr;\n"
            f"{fullInd}{ind}}}\n"
            f"{fullInd}}}\n"
        )

        # Write toYaml function
        functionname = f"{self.namespace}::{self.classname}::toYaml"
        target.write(
            f"{fullInd}inline auto {functionname}([[maybe_unused]] "
            f"::{common_namespace}::store_config const& config) const -> YAML::Node {{\n"
            f"{fullInd}{ind}using ::{common_namespace}::toYaml;\n"
            f"{fullInd}{ind}return toYaml(*value, config);\n"
            f"{fullInd}}}\n"
        )

        # Write fromYaml function
        functionname = f"{self.namespace}::{self.classname}::fromYaml"
        target.write(
            f"{fullInd}inline void {functionname}([[maybe_unused]] YAML::Node const& n) {{\n"
            f"{fullInd}{ind}using ::{common_namespace}::fromYaml;\n"
            f"{fullInd}{ind}fromYaml(n, *value);\n"
            f"{fullInd}}}\n"
        )


class EnumDefinition:
    """Prototype of a enum."""

    def __init__(self, name: str, values: list[str]):
        """Initialize enum definition with a name and possible values."""
        self.name = name
        self.values = values

        (self.namespace, self.classname) = split_name(name)
        self.namespace = safenamespacename(self.namespace)
        self.classname = safename(self.classname)

    def writeDefinition(self, target: IO[str], ind: str, common_namespace: str) -> None:
        """Write enum definition to output."""
        namespace = ""
        if len(self.name.split("#")) == 2:
            (namespace, classname) = split_name(self.name)
            namespace = safenamespacename(namespace)
            classname = safename(classname)

            name = namespace + "::" + classname
        else:
            name = safename(self.name)
            classname = name
        if len(namespace) > 0:
            target.write(f"namespace {namespace} {{\n")
        target.write(f"enum class {classname} : unsigned int {{\n{ind}")
        target.write(f",\n{ind}".join(map(safename, self.values)))
        target.write("\n};\n")
        target.write(f"inline auto to_string({classname} v) {{\n")
        target.write(f"{ind}static auto m = std::vector<std::string_view> {{\n")
        target.write(f'{ind}    "')
        target.write(f'",\n{ind}    "'.join(self.values))
        target.write(f'"\n{ind}}};\n')

        target.write(f"{ind}using U = std::underlying_type_t<{name}>;\n")
        target.write(f"{ind}return m.at(static_cast<U>(v));\n}}\n")

        if len(namespace) > 0:
            target.write("}\n")

        target.write(f"inline void to_enum(std::string_view v, {name}& out) {{\n")
        target.write(f"{ind}static auto m = std::map<std::string, {name}, std::less<>> {{\n")
        for v in self.values:
            target.write(f"""{ind}{ind}{{{q(v)}, {name}::{safename(v)}}},\n""")
        target.write(f"{ind}}};\n{ind}auto iter = m.find(v);\n")
        target.write(f"{ind}if (iter == m.end()) throw bool{{}};\n")
        target.write(f"{ind}out = iter->second;\n}}\n")

        # Write toYaml function
        target.write(f"namespace {common_namespace} {{\n")
        target.write(
            f"inline auto toYaml({name} v, [[maybe_unused]] "
            f"::{common_namespace}::store_config const& config) {{\n"
        )
        target.write(f"{ind}auto n = YAML::Node{{std::string{{to_string(v)}}}};\n")
        target.write(f'{ind}if (config.generateTags) n.SetTag("{name}");\n')
        target.write(f"{ind}return n;\n}}\n")

        # Write fromYaml function
        target.write(f"inline void fromYaml(YAML::Node n, {name}& out) {{\n")
        target.write(f"{ind}to_enum(n.as<std::string>(), out);\n}}\n")

        if len(self.values):
            target.write(f"template <> struct IsConstant<{name}> : std::true_type {{}};\n")

        target.write("}\n")

        target.write("\n")


# !TODO way to many functions, most of these shouldn't exists
def isPrimitiveType(v: Any) -> bool:
    """Check if v is a primitive type."""
    if not isinstance(v, str):
        return False
    return v in ["null", "boolean", "int", "long", "float", "double", "string"]


def hasFieldValue(e: Any, f: str, v: Any) -> bool:
    """Check if e has a field f value."""
    if not isinstance(e, dict):
        return False
    if f not in e:
        return False
    return bool(e[f] in [v, f"https://w3id.org/cwl/salad#{v}"])


def isRecordSchema(v: Any) -> bool:
    """Check if v is of type record schema."""
    return hasFieldValue(v, "type", "record")


def isEnumSchema(v: Any) -> bool:
    """Check if v is of type enum schema."""
    if not hasFieldValue(v, "type", "enum"):
        return False
    if "symbols" not in v:
        return False
    if not isinstance(v["symbols"], list):
        return False
    return True


def isArray(v: Any) -> bool:
    """Check if v is of type array."""
    if not isinstance(v, list):
        return False
    for i in v:
        if not pred(i):
            return False
    return True


def pred(i: Any) -> bool:
    """Check if v is any of the simple types."""
    return (
        isPrimitiveType(i)
        or isRecordSchema(i)
        or isEnumSchema(i)
        or isArraySchema(i)
        or isMapSchema(i)
        or isUnionSchema(i)
        or isinstance(i, str)
    )


def isArraySchema(v: Any) -> bool:
    """Check if v is of type array schema."""
    if not hasFieldValue(v, "type", "array"):
        return False
    if "items" not in v:
        return False
    if not isinstance(v["items"], list):
        return False

    for i in v["items"]:
        if not (pred(i) or isArray(i)):
            return False
    return True


def isMapSchema(v: Any) -> bool:
    """Check if v is of type map schema."""
    if not hasFieldValue(v, "type", "map"):
        return False
    if "values" not in v:
        return False
    if not isinstance(v["values"], list):
        return False

    for i in v["values"]:
        if not (pred(i) or isArray(i)):
            return False
    return True


def isUnionSchema(v: Any) -> bool:
    """Check if v is of type union schema."""
    return hasFieldValue(v, "type", "union")


class CppCodeGen(CodeGenBase):
    """Generation of C++ code for a given Schema Salad definition."""

    def __init__(
        self,
        base: str,
        target: IO[str],
        examples: Optional[str],
        package: str,
        copyright: Optional[str],
        spdx_copyright_text: Optional[list[str]],
        spdx_license_identifier: Optional[str],
    ) -> None:
        """Initialize the C++ code generator."""
        super().__init__()
        self.base_uri = base
        self.target = target
        self.examples = examples
        self.package = package
        self.copyright = copyright
        self.spdx_copyright_text = spdx_copyright_text
        self.spdx_license_identifier = spdx_license_identifier

        self.classDefinitions: dict[str, ClassDefinition] = {}
        self.enumDefinitions: dict[str, EnumDefinition] = {}
        self.mapDefinitions: dict[str, MapDefinition] = {}
        self.unionDefinitions: dict[str, UnionDefinition] = {}
        self.documentRootTypes: list[ClassDefinition] = []

    def convertTypeToCpp(self, type_declaration: Union[list[Any], dict[str, Any], str]) -> str:
        """Convert a Schema Salad type to a C++ type."""
        if not isinstance(type_declaration, list):
            return self.convertTypeToCpp([type_declaration])

        if len(type_declaration) == 1:
            if type_declaration[0] in ("null", "https://w3id.org/cwl/salad#null"):
                return "std::monostate"
            elif type_declaration[0] in (
                "string",
                "http://www.w3.org/2001/XMLSchema#string",
            ):
                return "std::string"
            elif type_declaration[0] in ("int", "http://www.w3.org/2001/XMLSchema#int"):
                return "int32_t"
            elif type_declaration[0] in (
                "long",
                "http://www.w3.org/2001/XMLSchema#long",
            ):
                return "int64_t"
            elif type_declaration[0] in (
                "float",
                "http://www.w3.org/2001/XMLSchema#float",
            ):
                return "float"
            elif type_declaration[0] in (
                "double",
                "http://www.w3.org/2001/XMLSchema#double",
            ):
                return "double"
            elif type_declaration[0] in (
                "boolean",
                "http://www.w3.org/2001/XMLSchema#boolean",
            ):
                return "bool"
            elif type_declaration[0] == "https://w3id.org/cwl/salad#Any":
                return "std::any"
            elif type_declaration[0] == "https://w3id.org/cwl/cwl#Expression":
                return "cwl_expression_string"
            elif type_declaration[0] in (
                "PrimitiveType",
                "https://w3id.org/cwl/salad#PrimitiveType",
            ):
                return "std::variant<bool, int32_t, int64_t, float, double, std::string>"
            elif isinstance(type_declaration[0], dict):
                if "type" in type_declaration[0] and type_declaration[0]["type"] in (
                    "enum",
                    "https://w3id.org/cwl/salad#enum",
                ):
                    name = type_declaration[0]["name"]
                    if name not in self.enumDefinitions:
                        self.enumDefinitions[name] = EnumDefinition(
                            type_declaration[0]["name"],
                            list(map(shortname, type_declaration[0]["symbols"])),
                        )
                    if len(name.split("#")) != 2:
                        return safename(name)
                    (namespace, classname) = name.split("#")
                    return safenamespacename(namespace) + "::" + safename(classname)
                elif "type" in type_declaration[0] and type_declaration[0]["type"] in (
                    "array",
                    "https://w3id.org/cwl/salad#array",
                ):
                    items = type_declaration[0]["items"]
                    if isinstance(items, list):
                        ts = [self.convertTypeToCpp(i) for i in items]
                        name = ", ".join(ts)
                        return f"std::vector<std::variant<{name}>>"
                    else:
                        i = self.convertTypeToCpp(items)
                        return f"std::vector<{i}>"
                elif "type" in type_declaration[0] and type_declaration[0]["type"] in (
                    "map",
                    "https://w3id.org/cwl/salad#map",
                ):
                    values = type_declaration[0]["values"]
                    if isinstance(values, list):
                        ts = [self.convertTypeToCpp(i) for i in values]
                        name = ", ".join(ts)
                        return f"std::map<std::string, std::variant<{name}>>"
                    else:
                        i = self.convertTypeToCpp(values)
                        return f"std::map<std::string, {i}>"
                elif "type" in type_declaration[0] and type_declaration[0]["type"] in (
                    "record",
                    "https://w3id.org/cwl/salad#record",
                ):
                    n = type_declaration[0]["name"]
                    (namespace, classname) = split_name(n)
                    return safenamespacename(namespace) + "::" + safename(classname)

                n = type_declaration[0]["type"]
                (namespace, classname) = split_name(n)
                return safenamespacename(namespace) + "::" + safename(classname)

            if len(type_declaration[0].split("#")) != 2:
                _logger.debug(f"// something weird2 about {type_declaration[0]}")
                return cast(str, type_declaration[0])

            (namespace, classname) = split_name(type_declaration[0])
            return safenamespacename(namespace) + "::" + safename(classname)

        type_declaration = list(map(self.convertTypeToCpp, type_declaration))
        type_declaration = ", ".join(type_declaration)
        return f"std::variant<{type_declaration}>"

    def epilogue(self, root_loader: Optional[TypeDef]) -> None:
        """Trigger to generate the epilouge code."""
        # find common namespace

        common_namespace = os.path.commonprefix(
            list(
                map(
                    lambda x: x.namespace,
                    list(self.classDefinitions.values()) + list(self.enumDefinitions.values()),
                )
            )
        )
        common_namespace = re.sub("(::)+$", "", common_namespace)

        """Generate final part of our cpp file."""
        if self.spdx_copyright_text:
            for text in self.spdx_copyright_text:
                self.target.write(f"""// SPDX-FileCopyrightText: {text}\n""")

        if self.spdx_license_identifier:
            self.target.write(f"""// SPDX-License-Identifier: {self.spdx_license_identifier}\n""")
        self.target.write("#pragma once\n\n")

        self.target.write(
            """/* This file was generated using schema-salad code generator.
 *
 * The embedded document is subject to the license of the original schema.
 """
        )

        if self.copyright:
            self.target.write("* The original schema is {self.copyright}.\n")

        self.target.write("*/\n\n")

        self.target.write(
            """#include <any>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <optional>
#include <string>
#include <string_view>
#include <variant>
#include <vector>
#include <yaml-cpp/yaml.h>

"""
        )

        self.target.write(f"namespace {common_namespace} {{\n")
        self.target.write(
            """
struct store_config {
    bool simplifyTypes = true;
    bool transformListsToMaps = true;
    bool generateTags = false;
};

inline auto simplifyType(YAML::Node type, store_config const& config) -> YAML::Node {
    if (!config.simplifyTypes) return type;
    auto is_optional = [](YAML::Node const & node) {
        return node.IsSequence() && node.size() == 2u && node[0].Scalar() == "null";
    };

    auto is_array = [](YAML::Node const & node) {
        return node.IsMap() && node["type"].Scalar() == "array" && node["items"].IsScalar();
    };

    // 1. Collapsing optional scalar types into one option
    if (is_optional(type) && type[1].IsScalar()) {
        type = type[1].as<std::string>() + "?";
    }

    // 2. Collapsing array types into one option
    if (is_array(type)) {
        type = type["items"].as<std::string>() + "[]";
    }

    // 3. Collapsing optional array types into one option
    if (is_optional(type) && is_array(type[1])) {
        type = type[1]["items"].as<std::string>() + "[]?";
    }

    return type;
}

inline auto expandType(YAML::Node type) -> YAML::Node {
    auto ends_with = [](std::string str, std::string suffix) {
        if (str.size() < suffix.size()) return false;
        auto str_suffix = str.substr(str.size()-suffix.size(), suffix.size());
        return str_suffix == suffix;
    };

    // 0. If not a scalar type, nothing to do
    if (!type.IsDefined() || !type.IsScalar()) {
        return type;
    }

    auto str = type.as<std::string>();
    // 1. Check if optional array type and expand
    if (ends_with(str, "[]?")) {
        auto result = YAML::Node{};
        result.push_back(YAML::Node{"null"});
        auto array = YAML::Node{};
        array["type"] = "array";
        array["items"] = expandType(YAML::Node(str.substr(0, str.size()-3)));
        result.push_back(array);
        return result;
    }

    // 2. Expand array
    if (ends_with(str, "[]")) {
        auto array = YAML::Node{};
        array["type"] = "array";
        array["items"] = expandType(YAML::Node(str.substr(0, str.size()-2)));
        return array;
    }

    // 3. Expand optional scalar type
    if (ends_with(str, "?")) {
        auto result = YAML::Node{};
        result.push_back(YAML::Node{"null"});
        result.push_back(expandType(YAML::Node(str.substr(0, str.size()-1))));
        return result;
    }
    return type;
}

inline auto mergeYaml(YAML::Node n1, YAML::Node n2) {
    for (auto const& e : n2) {
        n1[e.first.as<std::string>()] = e.second;
    }
    return n1;
}

// declaring toYaml
inline auto toYaml(bool v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(float v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(double v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(char v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int8_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint8_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int16_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint16_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int32_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint32_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(int64_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(uint64_t v, [[maybe_unused]] store_config const&) { return YAML::Node{v}; }
inline auto toYaml(std::monostate const&, [[maybe_unused]] store_config const&) {
    return YAML::Node(YAML::NodeType::Undefined);
}
inline auto toYaml(std::string const& v, [[maybe_unused]] store_config const&) {
    return YAML::Node{v};
}

template <typename T, typename ...Args>
auto anyToYaml_impl(std::any const& a, [[maybe_unused]] store_config const& config) {
    if (auto v = std::any_cast<T const>(&a)) {
        return toYaml(*v, config);
    }
    if constexpr (sizeof...(Args) > 0) {
        return anyToYaml_impl<Args...>(a, config);
    }
    return toYaml(std::monostate{}, config);
}

inline auto toYaml(std::any const& a, [[maybe_unused]] store_config const& config) {
    return anyToYaml_impl<bool,
                          float,
                          double,
                          char,
                          int8_t,
                          uint8_t,
                          int16_t,
                          uint16_t,
                          int32_t,
                          uint32_t,
                          int64_t,
                          uint64_t,
                          std::string>(a, config);
}

// declaring fromYaml
inline void fromYaml(YAML::Node const& n, bool& v) {
    v = n.as<bool>();
}
inline void fromYaml(YAML::Node const& n, float& v) {
    v = n.as<float>();
}
inline void fromYaml(YAML::Node const& n, double& v) {
    v = n.as<double>();
}
inline void fromYaml(YAML::Node const& n, int32_t& v) {
    v = n.as<int32_t>();
}
inline void fromYaml(YAML::Node const& n, int64_t& v) {
    v = n.as<int64_t>();
}
inline void fromYaml(YAML::Node const& n, std::string& v) {
    v = n.as<std::string>();
}
inline void fromYaml(YAML::Node const&, std::any&) {
}
inline void fromYaml(YAML::Node const&, std::monostate&) {
}

inline void addYamlField(YAML::Node& node, std::string const& key, YAML::Node value) {
    if (value.IsDefined()) {
        node[key] = value;
    }
}

inline auto convertListToMap(YAML::Node list, std::string const& mapSubject,
                             std::string const& mapPredicate, store_config const& config) {
    if (!config.transformListsToMaps) return list;
    if (mapSubject.empty()) return list;
    if (list.size() == 0) return list;
    auto map = YAML::Node{};
    for (YAML::Node n : list) {
        auto key = n[mapSubject].as<std::string>();
        if (mapPredicate.empty() || n[mapPredicate].IsMap() || n.size() > 2) {
            n.remove(mapSubject);
            map[key] = n;
        } else {
            map[key] = n[mapPredicate];
        }
    }
    return map;
}
inline auto convertMapToList(YAML::Node map, std::string const& mapSubject,
                             std::string const& mapPredicate) {
    if (mapSubject.empty()) return map;
    if (!map.IsDefined()) return map;
    if (!map.IsMap()) return map;
    auto list = YAML::Node{};
    for (auto n : map) {
        if (mapPredicate.empty() || n.second.IsMap()) {
            n.second[mapSubject] = n.first;
            list.push_back(n.second);
        } else {
            auto n2 = YAML::Node{};
            n2[mapSubject] = n.first;
            n2[mapPredicate] = n.second;
            list.push_back(n2);
        }
    }
    return list;
}

template <typename T> struct IsConstant : std::false_type {};

// fwd declaring toYaml
template <typename T>
auto toYaml(std::vector<T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node;
template <typename T>
auto toYaml(std::map<std::string, T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node;
template <typename T>
auto toYaml(T const& t, [[maybe_unused]] store_config const& config) -> YAML::Node;
template <typename ...Args>
auto toYaml(std::variant<Args...> const& t, [[maybe_unused]] store_config const& config) -> YAML::Node;

// fwd declaring fromYaml
template <typename T>
void fromYaml(YAML::Node const& n, std::vector<T>& v);
template <typename T>
void fromYaml(YAML::Node const& n, std::map<std::string, T>& v);
template <typename T>
void fromYaml(YAML::Node const& n, T& t);
template <typename ...Args>
void fromYaml(YAML::Node const& n, std::variant<Args...>& t);

template <typename T>
struct DetectAndExtractFromYaml {
    auto operator()(YAML::Node const&) const -> std::optional<T> {
        return std::nullopt;
    }
};

// special cwl expression string
struct cwl_expression_string {
    std::string s;

    auto toYaml([[maybe_unused]] store_config const& config) const {
        auto n = YAML::Node{s};
        if (config.generateTags) {
            n.SetTag("Expression");
        }
        return n;
    }
    void fromYaml(YAML::Node const& n) {
        s = n.as<std::string>();
    }
};


template <>
struct DetectAndExtractFromYaml<std::monostate> {
    auto operator()(YAML::Node const& n) const -> std::optional<std::monostate> {
        if (!n.IsDefined()) return std::monostate{};
        return std::nullopt;
    }
};

template <typename S>
struct DetectAndExtractFromYaml_implScalar {
    auto operator()(YAML::Node const& n) const -> std::optional<S> {
        try {
            if (n.IsScalar()) return n.as<S>();
        } catch(...) {}
        return std::nullopt;
    }
};

template <> struct DetectAndExtractFromYaml<bool>        : DetectAndExtractFromYaml_implScalar<bool>{};
template <> struct DetectAndExtractFromYaml<float>       : DetectAndExtractFromYaml_implScalar<float>{};
template <> struct DetectAndExtractFromYaml<double>      : DetectAndExtractFromYaml_implScalar<double>{};
template <> struct DetectAndExtractFromYaml<int32_t>     : DetectAndExtractFromYaml_implScalar<int32_t>{};
template <> struct DetectAndExtractFromYaml<int64_t>     : DetectAndExtractFromYaml_implScalar<int64_t>{};
template <> struct DetectAndExtractFromYaml<std::string> : DetectAndExtractFromYaml_implScalar<std::string>{};

template <typename T>
struct DetectAndExtractFromYaml<std::vector<T>> {
    auto operator()(YAML::Node const& n) const -> std::optional<std::vector<T>> {
        if (!n.IsDefined()) return std::nullopt;
        if (!n.IsSequence()) return std::nullopt;
        auto res = std::vector<T>{};
        fromYaml(n, res);
        return res;
    }
};

template <typename T>
struct DetectAndExtractFromYaml<std::map<std::string, T>> {
    auto operator()(YAML::Node const& n) const -> std::optional<std::map<std::string, T>> {
        if (!n.IsDefined()) return std::nullopt;
        if (!n.IsMap()) return std::nullopt;
        auto res = std::map<std::string, T>{};
        fromYaml(n, res);
        return res;
    }
};

template <typename T>
class heap_object {
    std::unique_ptr<T> data = std::make_unique<T>();

public:
    using value_t = T;
    heap_object() noexcept(false) = default;
    heap_object(heap_object const& oth) {
        *data = *oth;
    }
    heap_object(heap_object&& oth) noexcept(noexcept(*data = std::move(*oth))) {
        *data = std::move(*oth);
    }

    template <typename T2>
    heap_object(T2 const& oth) {
        *data = oth;
    }
    template <typename T2>
    heap_object(T2&& oth) noexcept(noexcept(*data = std::forward<T2>(oth))) {
        *data = std::forward<T2>(oth);
    }

    ~heap_object();

    auto operator=(heap_object const& oth) -> heap_object& {
        *data = *oth;
        return *this;
    }
    auto operator=(heap_object&& oth) noexcept(noexcept(*data = std::move(*oth))) -> heap_object& {
        *data = std::move(*oth);
        return *this;
    }

    template <typename T2>
    auto operator=(T2 const& oth) -> heap_object& {
        *data = oth;
        return *this;
    }
    template <typename T2>
    auto operator=(T2&& oth) noexcept(noexcept(*data = std::forward<T2>(oth))) -> heap_object& {
        *data = std::forward<T2>(oth);
        return *this;
    }

    auto operator->() noexcept(true) -> T* {
        return data.get();
    }
    auto operator->() const noexcept(true) -> T const* {
        return data.get();
    }
    auto operator*() noexcept(true) -> T& {
        return *data;
    }
    auto operator*() const noexcept(true) -> T const& {
        return *data;
    }
};

}
"""
        )
        # main body, printing fwd declaration, class definitions, and then implementations

        for key in self.classDefinitions:
            self.classDefinitions[key].writeFwdDeclaration(self.target, "", "    ")
        for key in self.mapDefinitions:
            self.mapDefinitions[key].writeFwdDeclaration(self.target, "", "    ")
        for key in self.unionDefinitions:
            self.unionDefinitions[key].writeFwdDeclaration(self.target, "", "    ")

        # remove parent classes, that are specialized/templated versions
        for key in self.classDefinitions:
            if len(self.classDefinitions[key].specializationTypes) > 0:
                self.classDefinitions[key].extends = []

        # remove fields that are available in a parent class
        for key in self.classDefinitions:
            for field in self.classDefinitions[key].allfields:
                found = False
                for parent_key in self.classDefinitions[key].extends:
                    fullKey = parent_key["namespace"] + "#" + parent_key["classname"]
                    for f in self.classDefinitions[fullKey].allfields:
                        if f.name == field.name:
                            found = True
                            break
                    if found:
                        break

                if not found:
                    self.classDefinitions[key].fields.append(field)  # noqa: B038

        # write definitions
        for key in self.enumDefinitions:
            self.enumDefinitions[key].writeDefinition(self.target, "    ", common_namespace)
        for key in self.classDefinitions:
            self.classDefinitions[key].writeDefinition(self.target, "", "    ", common_namespace)
        for key in self.mapDefinitions:
            self.mapDefinitions[key].writeDefinition(self.target, "    ", common_namespace)
        for key in self.unionDefinitions:
            self.unionDefinitions[key].writeDefinition(self.target, "    ", common_namespace)

        # CPP23: std::unique_ptr in heap_object is constexpr.
        # Hence, the compiler will try to instantiate the destructor on definition.
        # If the destructor was defined inside heap_object, other classes would only
        # be forward declared at this point.
        # This results in an error, because the destructor cannot be generated for
        # incomplete types.
        # Therefore, the destructor is defined here, after all classes have been defined.
        self.target.write(
            f"namespace {common_namespace} {{\n"
            f"template <typename T> heap_object<T>::~heap_object() = default;\n}}\n\n"
        )

        # write implementations
        for key in self.classDefinitions:
            self.classDefinitions[key].writeImplDefinition(
                self.target, "", "    ", common_namespace
            )
        for key in self.mapDefinitions:
            self.mapDefinitions[key].writeImplDefinition(self.target, "", "    ", common_namespace)
        for key in self.unionDefinitions:
            self.unionDefinitions[key].writeImplDefinition(
                self.target, "", "    ", common_namespace
            )

        self.target.write(f"namespace {common_namespace} {{\n")
        self.target.write(
            """
template <typename T>
auto toYaml(std::vector<T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node {
    auto n = YAML::Node(YAML::NodeType::Sequence);
    for (auto const& e : v) {
        n.push_back(toYaml(e, config));
    }
    return n;
}

template <typename T>
auto toYaml(std::map<std::string, T> const& v, [[maybe_unused]] store_config const& config) -> YAML::Node {
    auto n = YAML::Node(YAML::NodeType::Map);
    for (auto const& [key, value] : v) {
        n[key] = toYaml(value, config);
    }
    return n;
}

template <typename T>
auto toYaml(T const& t, [[maybe_unused]] store_config const& config) -> YAML::Node {
    if constexpr (std::is_enum_v<T>) {
        return toYaml(t, config);
    } else {
        return t.toYaml(config);
    }
}

template <typename ...Args>
auto toYaml(std::variant<Args...> const& t, store_config const& config) -> YAML::Node {
    return std::visit([config](auto const& e) {
        return toYaml(e, config);
    }, t);
}

template <typename T>
void fromYaml(YAML::Node const& n, std::vector<T>& v){
    if (!n.IsSequence()) return;
    for (auto e : n) {
        v.emplace_back();
        fromYaml(e, v.back());
    }
}

template <typename T>
void fromYaml(YAML::Node const& n, std::map<std::string, T>& v){
    if (!n.IsMap()) return;
    for (auto e : n) {
        auto key = e.first.as<std::string>();
        fromYaml(e.second, v[key]);
    }
}

template <typename T>
void fromYaml(YAML::Node const& n, T& t){
    if constexpr (std::is_enum_v<T>) {
        fromYaml(n, t);
    } else {
        t.fromYaml(n);
    }
}

template <typename SomeVariant, typename Head, typename ...Args>
bool detectAndExtractFromYaml(YAML::Node const& n, SomeVariant& v, Head* = nullptr) {
    auto r = DetectAndExtractFromYaml<Head>{}(n);
    if (r) {
        v = *r;
        return true;
    }
    if constexpr (sizeof...(Args) > 0) {
        return detectAndExtractFromYaml<SomeVariant, Args...>(n, v);
    }
    return false;
}

template <typename SomeVariant, typename Head, typename Tail>
bool detectAndExtractFromYaml(YAML::Node const& n, std::variant<std::monostate, Tail>& v, Head* = nullptr) {
    auto r = DetectAndExtractFromYaml<Head>{}(n);
    if (r) {
        v = *r;
        return true;
    }
    auto t = Tail{};
    fromYaml(n, t);
    v = t;
    return true;
}

template <typename ...Args>
void fromYaml(YAML::Node const& n, std::variant<Args...>& v){
    bool found = detectAndExtractFromYaml<std::variant<Args...>, Args...>(n, v);
    if (!found) throw std::runtime_error{"didn't find any overload"};
}
"""
        )
        rootTypes = []
        for cd in self.documentRootTypes:
            rootTypes.append(f"{cd.namespace}::{cd.classname}")
        documentRootType = ", ".join(rootTypes)

        self.target.write(f"using DocumentRootType = std::variant<{documentRootType}>;")
        self.target.write(
            """
auto load_document_from_yaml(YAML::Node n) -> DocumentRootType {
    DocumentRootType root;
    fromYaml(n, root);
    return root;
}
auto load_document_from_string(std::string document) -> DocumentRootType {
    return load_document_from_yaml(YAML::Load(document));
}
auto load_document(std::filesystem::path path) -> DocumentRootType {
    return load_document_from_yaml(YAML::LoadFile(path.string()));
}
void store_document(DocumentRootType const& root, std::ostream& ostream, store_config config={}) {
    auto y = toYaml(root, config);

    YAML::Emitter out;
    out << y;
    ostream << out.c_str() << std::endl;
}
void store_document(DocumentRootType const& root, std::filesystem::path const& path, store_config config={}) {
    auto ofs = std::ofstream{path};
    store_document(root, ofs, config);
}
auto store_document_as_string(DocumentRootType const& root, store_config config={}) -> std::string {
    auto ss = std::stringstream{};
    store_document(root, ss, config);
    return ss.str();
}

}"""
        )

    def parseRecordField(self, field: dict[str, Any]) -> FieldDefinition:
        """Parse a record field."""
        (namespace, classname, fieldname) = split_field(field["name"])
        mapSubject = ""
        mapPredicate = ""
        typeDSL = False
        if "jsonldPredicate" in field:
            if "mapSubject" in field["jsonldPredicate"]:
                mapSubject = field["jsonldPredicate"]["mapSubject"]
            if "mapPredicate" in field["jsonldPredicate"]:
                mapPredicate = field["jsonldPredicate"]["mapPredicate"]
            if "typeDSL" in field["jsonldPredicate"]:
                typeDSL = field["jsonldPredicate"]["typeDSL"]

        if isinstance(field["type"], dict):
            if field["type"]["type"] == "enum":
                fieldtype = "Enum"
            else:
                fieldtype = self.convertTypeToCpp(field["type"])

        else:
            fieldtype = self.convertTypeToCpp(field["type"])

        return FieldDefinition(
            name=fieldname,
            typeStr=fieldtype,
            optional=False,
            mapSubject=mapSubject,
            mapPredicate=mapPredicate,
            typeDSL=typeDSL,
        )

    def parseRecordSchema(self, stype: dict[str, Any]) -> None:
        """Parse a record schema."""
        cd = ClassDefinition(name=stype["name"])
        cd.abstract = stype.get("abstract", False)

        if "extends" in stype:
            for ex in aslist(stype["extends"]):
                (base_namespace, base_classname) = split_name(ex)
                ext = {"namespace": base_namespace, "classname": base_classname}
                cd.extends.append(ext)

        if "specialize" in stype:
            for e in aslist(stype["specialize"]):
                cd.specializationTypes.append(e["specializeFrom"])

        if "fields" in stype:
            for field in stype["fields"]:
                cd.allfields.append(self.parseRecordField(field))

        self.classDefinitions[stype["name"]] = cd

        if stype.get("documentRoot", False):
            self.documentRootTypes.append(cd)

    def parseMapSchema(self, stype: dict[str, Any]) -> str:
        """Parse a map schema."""
        name = cast(str, stype["name"])
        if name not in self.mapDefinitions:
            self.mapDefinitions[name] = MapDefinition(
                name, list(map(self.convertTypeToCpp, stype["values"]))
            )
        return name

    def parseUnionSchema(self, stype: dict[str, Any]) -> str:
        """Parse a union schema."""
        name = cast(str, stype["name"])
        if name not in self.unionDefinitions:
            self.unionDefinitions[name] = UnionDefinition(
                name, list(map(self.convertTypeToCpp, stype["names"]))
            )
        return name

    def parseEnum(self, stype: dict[str, Any]) -> str:
        """Parse a schema salad enum."""
        name = cast(str, stype["name"])
        if name not in self.enumDefinitions:
            self.enumDefinitions[name] = EnumDefinition(
                name, list(map(shortname, stype["symbols"]))
            )
        return name

    def parse(self, items: list[dict[str, Any]]) -> None:
        """Parse sechema salad items.

        This function is being called from the outside and drives
        the whole code generation.
        """
        for stype in items:
            if "type" in stype and stype["type"] == "documentation":
                continue

            if not (pred(stype) or isArray(stype)):
                raise SchemaException("not a valid SaladRecordField")

            # parsing a record
            if isRecordSchema(stype):
                self.parseRecordSchema(stype)
            elif isMapSchema(stype):
                self.parseMapSchema(stype)
            elif isUnionSchema(stype):
                self.parseUnionSchema(stype)
            elif isEnumSchema(stype):
                self.parseEnum(stype)
            else:
                _logger.error(f"not parsed{stype}")

        self.epilogue(None)
        self.target.close()
