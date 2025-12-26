"""A schema language for describing JSON or YAML structured linked data documents."""

import logging
from typing import Final

__author__: Final = "peter.amstutz@curoverse.com"

_logger: Final = logging.getLogger("salad")
_logger.addHandler(logging.StreamHandler())
_logger.setLevel(logging.INFO)
