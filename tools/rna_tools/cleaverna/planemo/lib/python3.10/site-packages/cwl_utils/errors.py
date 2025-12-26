# SPDX-License-Identifier: Apache-2.0
"""Common Exceptions."""


class ArrayMissingItems(BaseException):
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""


class JavascriptException(Exception):
    pass


class MissingKeyField(BaseException):
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""


class MissingTypeName(BaseException):
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""


class RecordMissingFields(BaseException):
    """From https://github.com/rabix/sbpack/blob/b8404a0859ffcbe1edae6d8f934e51847b003320/sbpack/lib.py ."""


class SubstitutionError(Exception):
    pass


class WorkflowException(Exception):
    pass


class GraphTargetMissingException(WorkflowException):
    """When a $graph is encountered and there is no target and no main/#main."""
