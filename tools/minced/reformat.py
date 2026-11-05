#!/usr/bin/env python

"""
    Extract importent information from the standard output file and put it in some standard format, like BED and tabular.
"""

import sys

bed = open(sys.argv[2], "w+")
tabular = open(sys.argv[3], "w+")

for line in open(sys.argv[1]):
    # Sequence 'CRISPRs' (10798 bp)
    if line.startswith("Sequence "):
        organism = line.split("'")[1]
    # CRISPR 1   Range: 679197 - 682529
    if line.startswith("CRISPR "):
        start, end = line.split("Range:")[1].strip().split("-")
        start = start.strip()
        end = end.strip()
        # minced reports 1-based inclusive coordinates, while BED uses zero-based
        # half-open intervals, so the start is shifted by one and the end is kept.
        bed.write("%s\t%d\t%s\n" % (organism, int(start) - 1, end))
    if line.rstrip().endswith("]"):
        cols = line.split()
        tabular.write(
            "%s\t%s\t%s\t%s\t%s\t%s\n"
            % (organism, cols[0], cols[1], cols[2], cols[4].rstrip(","), cols[5])
        )
