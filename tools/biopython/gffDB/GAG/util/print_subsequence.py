#!/usr/bin/env python

import sys

if len(sys.argv) != 5:
    sys.stderr.write("usage: python print_subsequence.py <fasta file> <sequence id> <start base> <end base>\n")
    sys.exit()

fasta_file = sys.argv[1]
seq_id = sys.argv[2]
start = int(sys.argv[3])
end = int(sys.argv[4])

with open(fasta_file, 'r') as fasta:
    found_seq_id = False
    sequence = ""
    for line in fasta:
        if line.startswith(">"):
            if seq_id == line.strip()[1:]:
                found_seq_id = True
                sys.stdout.write(">" + seq_id + " subsequence from " + str(start) + " to " + str(end) + "\n")
        else:
            if found_seq_id:
                sequence += line.strip()

sys.stdout.write(sequence[start-1:end] + "\n")        
