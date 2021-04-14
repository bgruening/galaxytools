#!/usr/bin/env python
"""
    Calculate the frequency of the nucleotides, the length and the GC-Content of a fastafile
    Bjoern Gruening
"""


import sys, os
from Bio import SeqIO


def __main__():
    input_file = open(sys.argv[1], "r")
    output_file = open(sys.argv[2], "w")

    output_file.write("Gene\tA\tC\tG\tT\tLength\tCG%\n")

    for record in SeqIO.parse(input_file, "fasta"):
        gene_name = record.name

        Acount = record.seq.count("A") + record.seq.count("a")
        Ccount = record.seq.count("C") + record.seq.count("c")
        Gcount = record.seq.count("G") + record.seq.count("g")
        Tcount = record.seq.count("T") + record.seq.count("t")

        seq_length = len(record.seq)

        GCcontent = float(Ccount + Gcount) / seq_length

        output_line = "%s\t%i\t%i\t%i\t%i\t%i\t%f\n" % (
            gene_name,
            Acount,
            Ccount,
            Gcount,
            Tcount,
            seq_length,
            GCcontent,
        )

        output_file.write(output_line)

    output_file.close()
    input_file.close()


if __name__ == "__main__":
    __main__()
