#! /usr/bin/env python

# Does its very best to parse an NCBI genome submission error report
# containing locations of regions to trim
# and write it to a bed file

import sys

if len(sys.argv) != 2:
    sys.stderr.write("usage: python error_report_to_bed.py <report.txt>\n")
    sys.exit()

report_file = sys.argv[1]
regions = {}

def parse_regions(text):
    """Return a list of (start, end) tuples."""
    regions = []
    region_pairs = text.strip().split(",")
    for region_pair in region_pairs:
        split_pair = region_pair.split("..")
        start = split_pair[0]
        end = split_pair[1]
        regions.append( [start, end] )
    return regions

with open(report_file, 'r') as report:
    for line in report:
        fields = line.strip().split("\t")
        if len(fields) != 4:
            # Hopefully it's a comment or something
            continue
        if "BioProject" in line or "PRJNA" in line:
            # Header line
            continue
        seq = fields[0].strip()
        regions = parse_regions(fields[2])
        for pair in regions:
            sys.stdout.write("\t".join( [seq, pair[0], pair[1]] ) + "\n")
