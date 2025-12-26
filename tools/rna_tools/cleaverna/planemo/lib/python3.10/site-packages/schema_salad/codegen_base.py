"""Base class for the generation of loaders from schema-salad definitions."""

from collections import OrderedDict
from collections.abc import MutableSequence
from typing import Any, Final, Optional, Union


class TypeDef:  # pylint: disable=too-few-public-methods
    """Schema Salad type description."""

    __slots__: Final = [
        "name",
        "init",
        "is_uri",
        "scoped_id",
        "ref_scope",
        "loader_type",
        "instance_type",
        "abstract",
    ]

    # switch to class-style typing.NamedTuple once support for Python < 3.6
    # is dropped
    def __init__(
        self,  # pylint: disable=too-many-arguments
        name: str,
        init: str,
        is_uri: bool = False,
        scoped_id: bool = False,
        ref_scope: Optional[int] = 0,
        loader_type: Optional[str] = None,
        instance_type: Optional[str] = None,
        abstract: bool = False,
    ) -> None:
        self.name: Final = name
        self.init: Final = init
        self.is_uri: Final = is_uri
        self.scoped_id: Final = scoped_id
        self.ref_scope: Final = ref_scope
        self.abstract: Final = abstract
        # Follow attributes used by Java but not Python.
        self.loader_type: Final = loader_type
        self.instance_type: Final = instance_type


class LazyInitDef:
    """Lazy initialization logic."""

    __slots__: Final = (
        "name",
        "init",
    )

    def __init__(
        self,
        name: str,
        init: str,
    ) -> None:
        """Create a LazyInitDef object."""
        self.name: Final = name
        self.init: Final = init


class CodeGenBase:
    """Abstract base class for schema salad code generators."""

    def __init__(self) -> None:
        self.collected_types: OrderedDict[str, TypeDef] = OrderedDict()
        self.lazy_inits: OrderedDict[str, LazyInitDef] = OrderedDict()
        self.vocab: dict[str, str] = {}

    def declare_type(self, declared_type: TypeDef) -> TypeDef:
        """Add this type to our collection, if needed."""
        if declared_type not in self.collected_types.values():
            self.collected_types[declared_type.name] = declared_type
        return declared_type

    def add_lazy_init(self, lazy_init: LazyInitDef) -> None:
        """Add lazy initialization logic for a given type."""
        self.lazy_inits[lazy_init.name] = lazy_init

    def add_vocab(self, name: str, uri: str) -> None:
        """Add the given name as an abbreviation for the given URI."""
        self.vocab[name] = uri

    def prologue(self) -> None:
        """Trigger to generate the prolouge code."""
        raise NotImplementedError()

    @staticmethod
    def safe_name(name: str) -> str:
        """Generate a safe version of the given name."""
        raise NotImplementedError()

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
        raise NotImplementedError()

    def end_class(self, classname: str, field_names: list[str]) -> None:
        """Signal that we are done with this class."""
        raise NotImplementedError()

    def type_loader(
        self,
        type_declaration: Union[list[Any], dict[str, Any]],
        container: Optional[str] = None,
        no_link_check: Optional[bool] = None,
    ) -> TypeDef:
        """Parse the given type declaration and declare its components."""
        raise NotImplementedError()

    def declare_field(
        self,
        name: str,
        fieldtype: TypeDef,
        doc: Optional[str],
        optional: bool,
        subscope: Optional[str],
    ) -> None:
        """Output the code to load the given field."""
        raise NotImplementedError()

    def declare_id_field(
        self,
        name: str,
        fieldtype: TypeDef,
        doc: Optional[str],
        optional: bool,
    ) -> None:
        """Output the code to handle the given ID field."""
        raise NotImplementedError()

    def uri_loader(
        self,
        inner: TypeDef,
        scoped_id: bool,
        vocab_term: bool,
        ref_scope: Optional[int],
        no_link_check: Optional[bool] = None,
    ) -> TypeDef:
        """Construct the TypeDef for the given URI loader."""
        raise NotImplementedError()

    def idmap_loader(
        self, field: str, inner: TypeDef, map_subject: str, map_predicate: Optional[str]
    ) -> TypeDef:
        """Construct the TypeDef for the given mapped ID loader."""
        raise NotImplementedError()

    def typedsl_loader(self, inner: TypeDef, ref_scope: Optional[int]) -> TypeDef:
        """Construct the TypeDef for the given DSL loader."""
        raise NotImplementedError()

    def epilogue(self, root_loader: TypeDef) -> None:
        """Trigger to generate the epilouge code."""
        raise NotImplementedError()

    def secondaryfilesdsl_loader(self, inner: TypeDef) -> TypeDef:
        """Construct the TypeDef for secondary files."""
        raise NotImplementedError()
