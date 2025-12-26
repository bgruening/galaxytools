# Copyright (C) 2017-2024 Vanessa Sochat.

# This Source Code Form is subject to the terms of the
# Mozilla Public License, v. 2.0. If a copy of the MPL was not distributed
# with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


def decodeUtf8String(inputStr):
    """Convert an UTF8 sequence into a string

    Required for compatibility with Python 2 where str==bytes
    inputStr -- Either a str or bytes instance with UTF8 encoding
    """
    return (
        inputStr
        if isinstance(inputStr, str) or not isinstance(inputStr, bytes)
        else inputStr.decode("utf8")
    )
