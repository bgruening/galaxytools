#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def process_unlabeled_file(input_file, label):
    with open(input_file, 'r') as handle:
        for record in SeqIO.parse(handle, 'fasta'):
            # Modify the description to include the label
            new_description = f"{record.description} LABEL={label}"
            yield SeqRecord(record.seq, id=record.id, description=new_description)

def main():
    parser = argparse.ArgumentParser(description="Preprocess and merge unlabeled FASTA files")
    parser.add_argument("--output", required=True, help="Output merged FASTA file")
    
    # Parse known args first
    args, unknown = parser.parse_known_args()
    
    # Process remaining arguments
    input_files = []
    for arg in unknown:
        if arg.startswith('--input_positive_') or arg.startswith('--input_negative_'):
            input_files.append((arg, next(iter(unknown[unknown.index(arg)+1:]), None)))

    # Process and write all sequences to the output file
    with open(args.output, 'w') as output_handle:
        for arg, file_path in input_files:
            if arg.startswith('--input_positive_'):
                label = 1
            elif arg.startswith('--input_negative_'):
                label = 0
            else:
                continue  # Skip any unrecognized arguments
            
            if file_path and not file_path.startswith('--'):
                SeqIO.write(process_unlabeled_file(file_path, label), output_handle, 'fasta')

if __name__ == "__main__":
    main()