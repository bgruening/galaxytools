"""Top level of GFF parsing providing shortcuts for useful classes.
"""
from GFFOutput import GFF3Writer, write  # noqa F401
from GFFParser import (DiscoGFFParser, GFFExaminer, GFFParser,  # noqa F401
                       parse, parse_simple)
