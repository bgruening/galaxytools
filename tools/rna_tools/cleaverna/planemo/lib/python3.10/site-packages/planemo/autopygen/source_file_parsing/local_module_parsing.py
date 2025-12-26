"""
Module responsible for resolving assignments and constant list comprehensives
used in argument parser
"""

import ast
import logging
from typing import (
    Any,
    List,
    Set,
    Tuple,
)

from planemo.autopygen.source_file_parsing.constants import WARNING_STRING
from planemo.autopygen.source_file_parsing.parsing_commons import add_parents
from planemo.autopygen.source_file_parsing.parsing_exceptions import CouldNotFixNameError


class UnknownNamesRemoval(ast.NodeVisitor):
    """
    Removes unknown names that can't be resolved and replaces them with a
    constant string that can be detected by linter

    Attributes
    ---

    unknown: Set[str]
      set of names that represent variables used in parser initialization,
      that have not been resolved yet

    """

    def __init__(self, unknown: Set[str]):
        self.unknown = unknown

    # currently able to resolve add_argument calls containing unknown names,
    # and list comprehension assignment
    def _reach_top(self, node: ast.Name) -> Tuple[ast.Call, ast.Name]:
        current = node
        parent = current.parent  # type: ignore

        def _reach_add_argument():
            return (
                isinstance(parent, ast.Call)
                and isinstance(parent.func, ast.Attribute)
                and parent.func.attr == "add_argument"
            )

        def _reach_assignment_as_list_comprehension():
            return isinstance(parent, ast.Assign) and isinstance(current, ast.ListComp)

        while not (_reach_add_argument() or _reach_assignment_as_list_comprehension()):
            current = parent
            if not hasattr(current, "parent"):
                raise CouldNotFixNameError
            parent = current.parent  # type: ignore

        return parent, current

    def _fix_name(self, node: ast.Name, name: str) -> bool:
        try:
            parent, current = self._reach_top(node)
        except CouldNotFixNameError:
            return False

        not_found_const = ast.Constant(value=f"{WARNING_STRING} Name {name}" f" could not be loaded")
        # if top is assignment
        if isinstance(parent, ast.Assign):
            logging.warning(f"Problem with assignment to {parent.targets[0].id}")  # type: ignore
            parent.value = ast.List(elts=[not_found_const], ctx=ast.Load())
            return True

        # this name can be a part of normal args
        if current in parent.args:
            idx = parent.args.index(current)  # type: ignore

            parent.args[idx] = not_found_const
            return True

        # or a part of keyword args
        current.value = not_found_const  # type: ignore
        return True

    # we dont care about class defifnitions, they should only depend
    # on outside things
    def visit_ClassDef(self, node: ast.ClassDef) -> Any:
        return

    def visit_Name(self, node: ast.Name) -> Any:
        if node.id not in self.unknown:
            return

        if not self._fix_name(node, node.id):
            logging.error(
                "Name could not be fixed and it can't be"
                " replaced with a constant,"
                " the parser extraction process failed"
            )


def handle_local_module_names(actions: List[ast.stmt], unknown_names: Set[str]) -> ast.Module:
    """
    Function used to remove assignments and list comprehensions which can't be
    resolved

    Parameters
    ----------
    actions: List[ast.stmt]
      list of actions extracted so far
    unknown_names:
      set of unknown names that have to be extracted

    Returns
    -------
    Python module containing assignments and list comprehensions whose values
    are based on constants.py
    """
    module = ast.Module(body=actions, type_ignores=[])
    add_parents(module)

    removal = UnknownNamesRemoval(unknown_names)
    removal.visit(module)

    return module
