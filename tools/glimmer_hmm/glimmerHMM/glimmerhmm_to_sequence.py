#!/usr/bin/env python
"""Convert GlimmerHMM gene predictions into protein sequences.

This works with both the GFF and the costumn Tabular Output.
And is only a wrapper to call the appropiate scripts.

Usage:
    glimmerhmm_to_sequence.py <glimmer output> <ref fasta> <output file> <format> <protein>

"""
import os
import sys

import glimmerhmm_gff_to_sequence
import glimmerhmm_tabular_to_sequence


def main(glimmer_file, ref_file, out_file, to_protein=False):
    if to_protein == "True":
        to_protein = True
    else:
        to_protein = False

    glimmerhmm_gff_to_sequence.main(glimmer_file, ref_file, out_file, to_protein)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print __doc__
        sys.exit()
    main(*sys.argv[1:])
