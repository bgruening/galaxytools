"""Support Singularity{,-CE} {2,3}.x or Apptainer 1.x."""

import re
from subprocess import check_output  # nosec
from typing import Optional

from .loghandler import _logger

# Cached version number of singularity
# This is a list containing major and minor versions as integer.
# (The number of minor version digits can vary among different distributions,
#  therefore we need a list here.)
_SINGULARITY_VERSION: Optional[list[int]] = None
# Cached flavor / distribution of singularity
# Can be singularity, singularity-ce or apptainer
_SINGULARITY_FLAVOR: str = ""


def get_version() -> tuple[list[int], str]:
    """
    Parse the output of 'singularity --version' to determine the flavor and version.

    Both pieces of information will be cached.

    :returns: A tuple containing:
              - A tuple with major and minor version numbers as integer.
              - A string with the name of the singularity flavor.
    """
    global _SINGULARITY_VERSION  # pylint: disable=global-statement
    global _SINGULARITY_FLAVOR  # pylint: disable=global-statement
    if _SINGULARITY_VERSION is None:
        version_output = check_output(  # nosec
            ["singularity", "--version"], text=True
        ).strip()

        version_match = re.match(r"(.+) version ([0-9\.]+)", version_output)
        if version_match is None:
            raise RuntimeError("Output of 'singularity --version' not recognized.")

        version_string = version_match.group(2)
        _SINGULARITY_VERSION = [int(i) for i in version_string.split(".")]
        _SINGULARITY_FLAVOR = version_match.group(1)

        _logger.debug(
            f"Singularity version: {version_string}" " ({_SINGULARITY_FLAVOR}."
        )
    return (_SINGULARITY_VERSION, _SINGULARITY_FLAVOR)


def is_apptainer_1_or_newer() -> bool:
    """
    Check if apptainer singularity distribution is version 1.0 or higher.

    Apptainer v1.0.0 is compatible with SingularityCE 3.9.5.
    See: https://github.com/apptainer/apptainer/releases
    """
    v = get_version()
    if v[1] != "apptainer":
        return False
    return v[0][0] >= 1


def is_version_2_6() -> bool:
    """
    Check if this singularity version is exactly version 2.6.

    Also returns False if the flavor is not singularity or singularity-ce.
    """
    v = get_version()
    if v[1] != "singularity" and v[1] != "singularity-ce":
        return False
    return v[0][0] == 2 and v[0][1] == 6


def is_version_3_or_newer() -> bool:
    """Check if this version is singularity version 3 or newer or equivalent."""
    if is_apptainer_1_or_newer():
        return True  # this is equivalent to singularity-ce > 3.9.5
    v = get_version()
    return v[0][0] >= 3


def is_version_3_1_or_newer() -> bool:
    """Check if this version is singularity version 3.1 or newer or equivalent."""
    if is_apptainer_1_or_newer():
        return True  # this is equivalent to singularity-ce > 3.9.5
    v = get_version()
    return v[0][0] >= 4 or (v[0][0] == 3 and v[0][1] >= 1)


def is_version_3_4_or_newer() -> bool:
    """Detect if Singularity v3.4+ is available."""
    if is_apptainer_1_or_newer():
        return True  # this is equivalent to singularity-ce > 3.9.5
    v = get_version()
    return v[0][0] >= 4 or (v[0][0] == 3 and v[0][1] >= 4)
