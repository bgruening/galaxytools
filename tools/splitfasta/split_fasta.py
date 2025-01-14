#!/usr/bin/env python

import argparse
import os

from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--records", type=int, default=None)
parser.add_argument("--limit", type=int, default=None)
parser.add_argument("--num-chunks", type=int, default=0)
parser.add_argument("input_file")
args = parser.parse_args()

input_filename = args.input_file
num_chunks = args.num_chunks
record_count = args.records
record_limit = args.limit

os.mkdir("splits")

if record_limit and num_chunks > record_limit:
    exit(f"ERROR: Requested number of chunks {num_chunks} exceeds limit {record_limit}")

if not record_count and (num_chunks != 0 or record_limit):
    # if no count is provided and if splitting into chunks or a limit is set, we need to count how many records are in the input file
    record_count = 0
    with open(input_filename) as input_file:
        for line in input_file:
            if line.lstrip().startswith(">"):
                record_count += 1

if num_chunks != 0:
    records_per_chunk = round(float(record_count) / num_chunks)

if record_limit and record_count > record_limit:
    exit(f"ERROR: Number of sequences {record_count} exceeds limit {record_limit}")

count = 1
with open(input_filename) as input_file:

    chunk_record_count = 0  # how many lines have we written to the output file
    records = []
    for record in SeqIO.parse(input_file, "fasta"):
        records.append(record)
        if num_chunks == 0 or (
            count < num_chunks and len(records) >= records_per_chunk
        ):
            if num_chunks == 0:
                output_filename = os.path.join("splits", record.id)
            else:
                output_filename = os.path.join("splits", "part{}".format(count))
            SeqIO.write(records, output_filename, "fasta")
            count += 1
            records = []

    if records:
        # this only applies for the mode where input file is
        # split into chunks
        output_filename = os.path.join("splits", "part{}".format(count))
        SeqIO.write(records, output_filename, "fasta")
