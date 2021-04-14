#!/usr/bin/env python
"""
Input: DNA Fasta File + Glimmer ORF File
Output: ORF-sequence as FASTA file
Author: Bjoern Gruening
"""
import sys, os
import Bio.SeqIO
from Bio.SeqRecord import SeqRecord as record


def __main__():
    if len(sys.argv) >= 4:
        glimmerfile = sys.argv[1]
        sequence = sys.argv[2]
        orf2seq = sys.argv[3]
    else:
        print "Missing input values."
        sys.exit()
    glimmer2sequence(sequence, glimmerfile, orf2seq)


def glimmer2sequence(
    sequence_path, glimmer_path, output_path, to_protein=False, translation_table=11
):

    fastafile = Bio.SeqIO.parse(open(sequence_path), "fasta")
    glimmerfile = open(glimmer_path, "r")
    orf2seq = open(output_path, "w")

    sequences = {}
    for entry in fastafile:
        sequences[entry.description] = entry

    for line in glimmerfile:
        if line.startswith(">"):
            entry = sequences[line[1:].strip()]
        else:
            columns = line.strip("\t").split()
            try:
                orf_start = int(columns[1])
                orf_end = int(columns[2])
            except:
                sys.stderr.write(
                    "Error: Failed to convert %s or %s to an integer. Is the input really a glimmer prediction file?\n"
                    % (columns[1], columns[2])
                )
                continue
            orf_name = columns[0]

            if orf_start <= orf_end:
                sequence = entry.seq[orf_start - 1 : orf_end]
                if to_protein:
                    sequence = sequence.translate(to_stop=True, table=translation_table)
                new_line = (
                    record(sequence, id=orf_name, description=entry.description).format(
                        "fasta"
                    )
                    + "\n"
                )
            else:
                sequence = entry.seq[orf_end - 1 : orf_start].reverse_complement()
                if to_protein:
                    sequence = sequence.translate(to_stop=True, table=translation_table)
                new_line = (
                    record(sequence, id=orf_name, description=entry.description).format(
                        "fasta"
                    )
                    + "\n"
                )
            orf2seq.write(new_line)

    orf2seq.close()
    glimmerfile.close()


if __name__ == "__main__":
    __main__()
