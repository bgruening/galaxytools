#!/usr/bin/env python
"""
Input: DNA Fasta File
Output: gezaehlte Codons ueber das gesamte File, nicht jedes einzelne ORF
Anika Erxleben
"""
import os
import sys

import Bio.SeqIO
import Bio.SeqUtils.CodonUsage as cu
from Bio.SeqRecord import SeqRecord as record

"""
strepto_orf2dna_out_plus1: anika/Annotation/data/strepto_orf2dna_out_plus1
codon_count: anika/Annotation/data/codon_count2  outputfile
"""


def __main__():

    if len(sys.argv) >= 3:
        orf2dna_out = sys.argv[1]
        codon_count = open(sys.argv[2], "w")

    else:
        print("da fehlt was da oben")
        sys.exit()

    a = cu.CodonAdaptationIndex()
    a.generate_index(orf2dna_out)

    # for record in Bio.SeqIO.parse(open(orf2dna_out), "fasta"):
    #    print record.seq

    for key, value in a.codon_count.items():
        codon_count.write(key + "\t" + str(value) + "\n")

    codon_count.close()


if __name__ == "__main__":
    __main__()
