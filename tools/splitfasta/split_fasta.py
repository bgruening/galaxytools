#!/usr/bin/env python

import os
import sys
from Bio import SeqIO

num_chunks = 0
if len(sys.argv) == 3:
    num_chunks = int(sys.argv[2])
    input_filename = sys.argv[1]
elif len(sys.argv) == 2:
    input_filename = sys.argv[1]
else:
    exit("Usage: split_fasta.py <input_filename> [<num_chunks>]")

os.mkdir('splits')

if num_chunks != 0:
    # if splitting into chunks we need to count how many records are in the
    # input file
    record_count = 0
    with open(input_filename) as input_file:
        for line in input_file:
            if line.lstrip().startswith('>'):
                record_count += 1

    records_per_chunk = round(float(record_count) / num_chunks)


count = 1
with open(input_filename) as input_file:

    chunk_record_count = 0  # how many lines have we written to the output file
    records = []
    for record in SeqIO.parse(input_file, 'fasta'):
        records.append(record)
        chunk_record_count += 1
        if num_chunks == 0 or (count < num_chunks and
           chunk_record_count >= records_per_chunk):
            if num_chunks == 0:
                output_filename = os.path.join('splits', record.id + '.fasta')
            else:
                output_filename = os.path.join('splits', 'part{}.fasta'.format(count))
            SeqIO.write(records, output_filename, 'fasta')
            count += 1
            records = []
            chunk_record_count = 0

    if records:
        # this only applies for the mode where input file is
        # split into chunks
        output_filename = os.path.join('splits', 'part{}.fasta'.format(count))
        SeqIO.write(records, output_filename, 'fasta')
