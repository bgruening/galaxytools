#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse
from Bio import SeqIO


def change_fasta_header(infile, outfile, mode):
    new_sequence_array = []

    ofile = open(outfile, 'w+')
    for line in open(infile):
        if line.startswith('>'):
            header = line[1:].strip()
            if mode == 'id_only':
                ofile.write('>%s\n' % header.split()[0])
            else:
                ofile.write('>%s\n' % header.replace(' ', '_'))
        else:
            ofile.write(line)

if __name__ == "__main__" :
    change_fasta_header(sys.argv[1], sys.argv[2], sys.argv[3])
