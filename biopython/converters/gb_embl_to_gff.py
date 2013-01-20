#!/usr/bin/env python
"""Convert a GenBank file into GFF format.

https://raw.github.com/chapmanb/bcbb/master/gff/Scripts/gff/genbank_to_gff.py

Usage:
    gb_embl_to_gff.py <genbank_file> <gff_file> <input_format embl or genbank>
"""
import sys
import os

from Bio import SeqIO
from Bio import Seq

from BCBio import GFF

def main(gb_file, ofile, iformat):
    with open(ofile, "w") as out_handle:
        GFF.write(SeqIO.parse(gb_file, iformat), out_handle)

if __name__ == "__main__":
    main(*sys.argv[1:])
