# SPDX-License-Identifier: Apache-2.0
"""Shared logging object."""
import logging

_logger = logging.getLogger("cwl_utils")  # pylint: disable=invalid-name
defaultStreamHandler = logging.StreamHandler()  # pylint: disable=invalid-name
_logger.addHandler(defaultStreamHandler)
_logger.setLevel(logging.INFO)
