#!/usr/bin/env python
import os
import sys
from Bio import SeqIO

if __name__ == "__main__":
    inpath = sys.argv[1]
    os.mkdir('splits')
    with open(inpath, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            header = os.path.join('splits', record.id + '.fasta')
            with open(header, 'w') as handle2:
                SeqIO.write([record], handle2, 'fasta')
