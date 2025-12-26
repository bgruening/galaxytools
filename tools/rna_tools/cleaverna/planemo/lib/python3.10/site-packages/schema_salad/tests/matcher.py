# Copyright (C) The Arvados Authors. All rights reserved.
#
# SPDX-License-Identifier: Apache-2.0

import difflib
import re
from typing import Any

from schema_salad.utils import json_dumps


class JsonDiffMatcher:
    """Raise AssertionError with a readable JSON diff when not __eq__().

    Used with assert_called_with() so it's possible for a human to see
    the differences between expected and actual call arguments that
    include non-trivial data structures.
    """

    def __init__(self, expected: Any):
        self.expected = expected

    def __eq__(self, actual: Any) -> bool:
        expected_json = json_dumps(self.expected, sort_keys=True, indent=2)
        actual_json = json_dumps(actual, sort_keys=True, indent=2)
        if expected_json != actual_json:
            raise AssertionError(
                "".join(
                    difflib.context_diff(
                        expected_json.splitlines(True),
                        actual_json.splitlines(True),
                        fromfile="Expected",
                        tofile="Actual",
                    )
                )
            )
        return True


def StripYAMLComments(yml: str) -> Any:
    return re.sub(r"(?ms)^(#.*?\n)*\n*", "", yml)
