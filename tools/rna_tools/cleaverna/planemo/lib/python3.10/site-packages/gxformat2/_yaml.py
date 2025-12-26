"""YAML loading utilities for gxformat2.

This was repeatedly used in Galaxy and Planemo so moved this module
to gxformat2.yaml.
"""
import warnings

from gxformat2.yaml import (
    ordered_dump,
    ordered_load,
)

__all__ = ('ordered_load', 'ordered_dump')

warnings.warn("Importing gxformat2._yaml is deprecated, use gxformat2.yaml instead", DeprecationWarning, stacklevel=2)
