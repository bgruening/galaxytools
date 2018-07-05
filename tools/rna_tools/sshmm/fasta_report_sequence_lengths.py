#!/usr/bin/python

import sys

"""
Input: FASTA file
Output: Print sequence ID and corresponding sequence length

Example output:
chr1	248956422
chr2	242193529
chr3	198295559
...

"""

# Check input.
if not len(sys.argv) == 2:
    exit("Usage: fasta_report_sequence_lengths.py <fasta_file>")

fasta_file = sys.argv[1]

seq_id = "id"
seq_len = 0

# Go through FASTA file, extract sequence lengths.
with open(fasta_file) as f:
    for line in f:
        if line.startswith(">"):
            new_id = line[1:].strip()
            if seq_len:
                print "%s\t%i" % (seq_id, seq_len)
            seq_len = 0
            seq_id = new_id
        else:
            seq_len += len(line.strip())
f.closed

# Print last sequence length.
if seq_len:
    print "%s\t%i" % (seq_id, seq_len)

