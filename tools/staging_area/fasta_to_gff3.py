#!/usr/bin/env python

import sys, os, re
from Bio import SeqIO

def main():

    input_file = open(sys.argv[1], "rU")
    output_file = open(sys.argv[2], 'w')

    filtered_records =[]
    output_file.write('##gff-version 3\n')
    for record in SeqIO.parse(input_file, "fasta") :
        seq_start = 1
        seq_id = record.description
        
        seq_end = seq_start + len(record.seq)
        output_file.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (
                (seq_id, 'scaffold_description', 'contig', seq_start, seq_end, '.', '.', '.', 'Name='+seq_id))
            )
        #seq_start = seq_end + 1
        #7180000016251	TBLASTN	contig	1	50000	.	.	.	Name=7180000016251


    SeqIO.write(filtered_records, output_file, "fasta")


if __name__ == "__main__":
    main()
