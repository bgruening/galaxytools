#!/usr/bin/env python

"""
    Extract importent information from the standard output file and put it in some standard format, like BED and tabular.
"""

import sys

bed = open(sys.argv[2], "w+")
tabular = open(sys.argv[3], "w+")

for line in open(sys.argv[1]):
    # ORGANISM:  gi|21226102|ref|NC_003901.1| Methanosarcina mazei Go1 chromosome, complete genome
    if line.startswith("ORGANISM:"):
        organism = line.lstrip("ORGANISM:").strip()
    # CRISPR 1   Range: 679197 - 682529
    if line.startswith("CRISPR "):
        start, end = line.split("Range:")[1].strip().split("-")
        start = start.strip()
        end = end.strip()
        bed.write("%s\t%s\t%s\n" % (organism, start, end))
    if line.rstrip().endswith("]"):
        cols = line.split()
        tabular.write(
            "%s\t%s\t%s\t%s\t%s\t%s\n"
            % (organism, cols[0], cols[1], cols[2], cols[4].rstrip(","), cols[5])
        )
