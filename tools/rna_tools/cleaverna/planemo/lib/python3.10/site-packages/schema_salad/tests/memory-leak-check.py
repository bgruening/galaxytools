#!/usr/bin/env python

import sys

import cwl_utils.parser  # type: ignore[import-not-found]
import objgraph  # type: ignore[import-untyped]

growth = []
for _ in range(5):
    growth = objgraph.growth(limit=50)
    cwl_utils.parser.load_document_by_uri(sys.argv[1])
objgraph.show_growth(limit=50)
if len(growth) != 0:
    exit(1)
