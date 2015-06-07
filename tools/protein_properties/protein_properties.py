#!/usr/bin/env python

import sys
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis

with open(sys.argv[2], 'w+') as result:

    result.write("ID\tSequence\tMW\tIP\tgravy\tlength\tinstability\tmonoisotpoic\n")
    for record in SeqIO.parse(open(sys.argv[1]), "fasta"):
        a = ProteinAnalysis(str(record.seq))

        properties = list()
        properties.append(record.id)
        properties.append(a.sequence)
        properties.append(a.molecular_weight())
        properties.append(a.isoelectric_point())
        properties.append(a.gravy())
        properties.append(a.length)
        properties.append(a.instability_index())
        properties.append(a.aromaticity())

        result.write( '\t'.join(map(str, properties))+"\n" )

