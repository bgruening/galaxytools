"""Shared Exception classes."""

from collections.abc import Sequence
from typing import Final, Optional, Union

from .sourceline import SourceLine, reflow_all, strip_duplicated_lineno


def _simplify(exc: "SchemaSaladException") -> list["SchemaSaladException"]:
    return [exc] if len(exc.message) else exc.children


def _with_bullet(exc: "SchemaSaladException", bullet: str) -> "SchemaSaladException":
    if exc.bullet == "":
        exc.bullet = bullet
    return exc


class SchemaSaladException(Exception):
    """Base class for all schema-salad exceptions."""

    def __init__(
        self,
        msg: str,
        sl: Optional[SourceLine] = None,
        children: Optional[Sequence["SchemaSaladException"]] = None,
        bullet_for_children: str = "",
        detailed_message: Optional[str] = None,
    ) -> None:
        super().__init__(msg)
        self.message: Final = self.args[0]
        self.detailed_message: Final = detailed_message
        self.file: Optional[str] = None
        self.start: Optional[tuple[int, int]] = None
        self.end: Optional[tuple[int, int]] = None

        self.is_warning: bool = False

        # It will be set by its parent
        self.bullet: str = ""

        if children is None:
            self.children: list["SchemaSaladException"] = []
        elif len(children) <= 1:
            self.children = sum((_simplify(c) for c in children), [])
        else:
            self.children = sum(
                (_simplify(_with_bullet(c, bullet_for_children)) for c in children), []
            )

        self.with_sourceline(sl)
        self.propagate_sourceline()

    def propagate_sourceline(self) -> None:
        if self.file is None:
            return
        for c in self.children:
            if c.file is None:
                c.file = self.file
                c.start = self.start
                c.end = self.end
                c.propagate_sourceline()

    def as_warning(self) -> "SchemaSaladException":
        self.is_warning = True
        for c in self.children:
            c.as_warning()
        return self

    def with_sourceline(self, sl: Optional[SourceLine]) -> "SchemaSaladException":
        """Use the provided SourceLine to set the causal location."""
        if sl and sl.file():
            self.file = sl.file()
            self.start = sl.start()
            self.end = sl.end()
        else:
            self.file = None
            self.start = None
            self.end = None
        return self

    def leaves(self) -> list["SchemaSaladException"]:
        """Return the list of all the exceptions at the tips of the tree."""
        if len(self.children) > 0:
            return sum((c.leaves() for c in self.children), [])
        if len(self.message):
            return [self]
        return []

    def prefix(self) -> str:
        pre: str = ""
        if self.file:
            linecol0: Union[int, str] = ""
            linecol1: Union[int, str] = ""
            if self.start:
                linecol0, linecol1 = self.start
            pre = f"{self.file}:{linecol0}:{linecol1}: "

        return pre + "Warning: " if self.is_warning else pre

    def summary(self, level: int = 0, with_bullet: bool = False) -> str:
        indent_per_level: Final = 2
        spaces: Final = (level * indent_per_level) * " "
        bullet: Final = self.bullet + " " if len(self.bullet) > 0 and with_bullet else ""
        message_string: Final = (
            self.detailed_message
            if (len(self.children) < 1 and self.detailed_message)
            else self.message
        )
        return f"{self.prefix()}{spaces}{bullet}{message_string}"

    def __str__(self) -> str:
        """Convert to a string using :py:meth:`pretty_str`."""
        return str(self.pretty_str())

    def pretty_str(self, level: int = 0) -> str:
        messages: Final = (
            len(self.message)
            if len(self.children) > 0
            else len(self.detailed_message or self.message)
        )
        my_summary: Final = [self.summary(level, True)] if messages else []
        next_level: Final = level + 1 if messages else level

        ret: Final = "\n".join(
            e for e in my_summary + [c.pretty_str(next_level) for c in self.children]
        )
        if level == 0:
            return strip_duplicated_lineno(reflow_all(ret))
        return ret


class SchemaException(SchemaSaladException):
    """Indicates error with the provided schema definition."""


class ValidationException(SchemaSaladException):
    """Indicates error with document against the provided schema."""


class ClassValidationException(ValidationException):
    pass


def to_one_line_messages(exc: SchemaSaladException) -> str:
    return "\n".join(c.summary() for c in exc.leaves())
