"""
Module responsible for discovery of import statements importing Argument parser
and discovery of the statements initializing the parser itself
"""

import ast
import logging
from typing import (
    Any,
    List,
    Optional,
    Set,
    Tuple,
)

from .constants import STD_LIB_MODULE_NAMES
from .decoy_parser import (
    CustomParserUnavailableException,
    obtain_class_def,
)
from .parsing_commons import Discovery
from .parsing_exceptions import (
    ArgParseImportNotFound,
    ArgParserNotUsed,
)

ARGPARSE_MODULE_NAME = "argparse"
ARGUMENT_PARSER_CLASS_NAME = "ArgumentParser"


def is_this_x_creation(node: ast.Assign, function_name: str):
    if not (len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)):
        return False, None

    name = node.targets[0].id
    if not (
        isinstance(node.value, ast.Call)
        and isinstance(node.value.func, ast.Attribute)
        and node.value.func.attr == function_name
    ):
        return False, None

    return True, name


class ImportDiscovery(Discovery):
    """
    Class responsible for discovery and extraction of import statements
    """

    def __init__(self, actions: List[ast.stmt]):
        super(ImportDiscovery, self).__init__(actions)
        self.argparse_module_alias: Optional[str] = None
        self.argument_parser_alias: Optional[str] = None
        self.known_names: Set[str] = set()

    def visit_Import(self, node: ast.Import) -> Any:
        for item in node.names:
            if item.name == ARGPARSE_MODULE_NAME:
                alias = item.asname if item.asname is not None else ARGPARSE_MODULE_NAME
                self.argparse_module_alias = alias
                self.known_names.add(alias)

            if item.name in STD_LIB_MODULE_NAMES:
                self.actions.append(node)

    def visit_ImportFrom(self, node: ast.ImportFrom) -> Any:
        if node.module is None:
            return

        for name in node.module.split("."):
            name_in_known_modules = name in STD_LIB_MODULE_NAMES
            if name_in_known_modules:
                self.actions.append(node)
                for item in node.names:
                    alias = item.asname or item.name
                    self.known_names.add(alias)
                    # in case argparse is being imported, determine the
                    # alias of the parser, if there is any
                    if name == ARGPARSE_MODULE_NAME and item.name == ARGUMENT_PARSER_CLASS_NAME:
                        self.argument_parser_alias = alias

    def report_findings(self) -> Tuple[List[ast.stmt], Optional[str], Optional[str], Set[str]]:
        if self.argparse_module_alias is None and self.argument_parser_alias is None:
            raise ArgParseImportNotFound("No argparse import found")

        return (self.actions, self.argparse_module_alias, self.argument_parser_alias, self.known_names)


class SimpleParserDiscoveryAndReplacement(Discovery):
    """
    Class responsible for discovery of ArgumentParser creation
    and assignment, and replacement of the class definition
    by the one supplied through constructor
    """

    def __init__(
        self, actions: List[ast.stmt], argparse_alias: str, argument_parser_alias: str, custom_parser_def: ast.ClassDef
    ):
        self.argument_parser_alias = argument_parser_alias
        self.argparse_module_alias = argparse_alias
        self.main_parser_name: Optional[str] = None
        self.argparse_found = False
        self.custom_parser_def = custom_parser_def

        super(SimpleParserDiscoveryAndReplacement, self).__init__(actions)

    @staticmethod
    def is_simple_assignment(node: ast.Assign):
        return len(node.targets) == 1 and isinstance(node.targets[0], ast.Name)

    @staticmethod
    def imported_using_from(node: ast.Assign, argument_parser_alias: str):
        return (
            isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Name)
            and node.value.func.id == argument_parser_alias
        )

    @staticmethod
    def imported_using_import(node: ast.Assign, argparse_module_alias: str):
        return (
            isinstance(node.value, ast.Call)
            and isinstance(node.value.func, ast.Attribute)
            and node.value.func.attr == ARGUMENT_PARSER_CLASS_NAME
            and node.value.func.value.id == argparse_module_alias  # type: ignore
        )

    def visit_Assign(self, node: ast.Assign):
        if self.argparse_found:
            return
        # visit into children of this node is not necessary
        if not self.is_simple_assignment(node):
            return

        name = node.targets[0].id  # type: ignore
        import_from = self.imported_using_from(node, self.argument_parser_alias)
        import_import = self.imported_using_import(node, self.argparse_module_alias)

        if import_import or import_from:
            self.actions.append(self.custom_parser_def)
            self.main_parser_name = name
            self.actions.append(node)
            self.argparse_found = True
            self._replace_parser(node, import_from)

    def _replace_parser(self, node: ast.Assign, imported_using_from: bool):
        # FIXME TODO currently, passing variables to custom argument parser
        # is not supported
        if node.value.args or node.value.keywords:  # type: ignore
            logging.warning(
                "Arguments that are normally passed to argument"
                " parser will be ignored. Their use is"
                " not currently supported"
            )
        node.value.args = []  # type: ignore
        node.value.keywords = []  # type: ignore
        if imported_using_from:
            self.custom_parser_def.bases[0] = ast.Name(self.argument_parser_alias, ast.Load())
            node.value.func.id = self.custom_parser_def.name  # type: ignore
            return

        assert type(node.value is ast.Call)
        node.value.func = ast.Name(self.custom_parser_def.name, ast.Load())  # type: ignore

    def report_findings(self) -> Tuple:
        if self.main_parser_name is None:
            raise ArgParserNotUsed("Argument parser not used")

        return self.actions, self.main_parser_name


GROUPS_AND_SUBPARSERS = ["add_argument_group", "add_parser", "add_subparsers"]


class GroupAndSubparsersDiscovery(Discovery):
    """
    Class responsible for discovery of statements that initialize argument
    groups and subparsers
    """

    def __init__(self, actions: List[ast.stmt], known_names: Set[str], main_name: str):
        self.main_name = main_name
        self.known_names = known_names
        super(GroupAndSubparsersDiscovery, self).__init__(actions)

    def visit_Assign(self, node: ast.Assign):
        for name in GROUPS_AND_SUBPARSERS:
            is_correct_creation, name = is_this_x_creation(node, name)
            if is_correct_creation:
                self.known_names.add(name)
                self.actions.append(node)

    def report_findings(self) -> Tuple:
        return self.actions, self.known_names


# # this visitor goes through all calls and extracts those to argument
# parser and groups. IMPORTANT! it also renames parsers on which those calls
# are called to ensure everything can be interpreted correctly
class ArgumentCreationDiscovery(Discovery):
    """
    Class responsible for extraction of statements which initialize the input
    arguments. It is able to extract function calls on the original parser,
    and on the argument groups extracted by GroupDiscovery
    """

    def __init__(self, actions: List[ast.stmt], main_name: str):
        self.main_name = main_name
        super(ArgumentCreationDiscovery, self).__init__(actions)

    @staticmethod
    def is_call_on_parser_or_group(node: ast.Call):
        return isinstance(node.func, ast.Attribute) and node.func.attr == "add_argument"

    def visit_Call(self, node: ast.Call) -> Any:
        if self.is_call_on_parser_or_group(node):
            self.actions.append(ast.Expr(node))

        self.generic_visit(node)

    def report_findings(self) -> Tuple[List[ast.stmt]]:
        return (self.actions,)


def get_parser_init_and_actions(source: ast.Module) -> Tuple[List[ast.stmt], str, Set[str]]:
    """
    Function used to extract necessary imports, parser and argument creation
     function calls

    Parameters
    ----------
    source : ast.Module
      source file parsed into ATT

    Returns
    -------
    List of extracted AST nodes, the main name of the parser and a set of
    section names
    """

    actions: List[ast.stmt] = []
    custom_parser_class_def = obtain_class_def()

    if custom_parser_class_def is None:
        raise CustomParserUnavailableException("Custom parser unavailable")

    import_discovery = ImportDiscovery(actions)
    actions, argparse_module_alias, argparse_class_alias, known_names = import_discovery.visit_and_report(source)
    known_names.add(custom_parser_class_def.name)

    parser_discovery = SimpleParserDiscoveryAndReplacement(
        actions, argparse_module_alias, argparse_class_alias, custom_parser_class_def
    )
    actions, parser_name = parser_discovery.visit_and_report(source)

    group_discovery = GroupAndSubparsersDiscovery(actions, known_names, parser_name)
    actions, known_names = group_discovery.visit_and_report(source)

    argument_creation = ArgumentCreationDiscovery(actions, parser_name)
    (actions,) = argument_creation.visit_and_report(source)

    return actions, parser_name, known_names
