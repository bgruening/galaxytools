#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

sys.stdout.write("ID\tMW\tIP\tgravy\tlength\tinstability\tmonoisotpoic\tSequence\n")

for record in SeqIO.parse(sys.stdin, "fasta"):
    a = ProteinAnalysis(str(record.seq))

    properties = list()
    properties.append(record.id)
    properties.append(a.molecular_weight())
    properties.append(a.isoelectric_point())
    properties.append(a.gravy())
    properties.append(a.length)
    properties.append(a.instability_index())
    properties.append(a.aromaticity())
    # always last column to make the output more readable
    properties.append(a.sequence)
    sys.stdout.write( '\t'.join(map(str, properties))+"\n" )

