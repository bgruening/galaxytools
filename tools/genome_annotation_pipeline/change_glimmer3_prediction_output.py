#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse


def change_glimmer3_prediction_output(infile, outfile):
    outfile = open(outfile, "w")
    contig_name = ""
    for line in open(infile):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith(">"):
            # Only take the Sequence ID
            contig_name = line[1:].split()[0]
            outfile.write(line + "\n")
        else:
            # orf-line, starts with the orfname at the beginning
            orf_name = line.split()[0]
            outfile.write(
                line.replace(orf_name, "%s_%s" % (contig_name, orf_name)) + "\n"
            )

    outfile.close()


if __name__ == "__main__":
    change_glimmer3_prediction_output(sys.argv[1], sys.argv[2])
