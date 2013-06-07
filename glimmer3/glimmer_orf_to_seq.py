#!/usr/bin/env python
"""
Input: DNA Fasta File + Glimmer ORF File
Output: ORF-sequence as FASTA file
Author: Bjoern Gruening
"""
import sys, os
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord as record

def __main__():
    if len(sys.argv) >= 4:
        glimmerfile = open(sys.argv [1], "r")
        sequence = open(sys.argv[2])
        orf2seq = open(sys.argv [3], "w")
    else:
        print "Missing input values."
        sys.exit()

    fastafile = Bio.SeqIO.parse(sequence, "fasta")

    sequences = {}
    for entry in fastafile:
        sequences[entry.description] = entry

    for line in glimmerfile:
        if line.startswith('>'):
            entry = sequences[ line[1:].strip() ]
        else:
            orf_start = int(line[8:17])
            orf_end = int(line[18:26])

            orf_name = line[0:8]
            if orf_start <= orf_end:
                new_line = record(entry.seq[orf_start-1 : orf_end], id = orf_name, description = entry.description).format("fasta") + "\n"
            else:         
                new_line = record(entry.seq[orf_end-1 : orf_start].reverse_complement(), id = orf_name, description = entry.description).format("fasta") + "\n"
            orf2seq.write(new_line)

    orf2seq.close()
    glimmerfile.close()
    sequence.close()

if __name__ == "__main__" :
    __main__()
