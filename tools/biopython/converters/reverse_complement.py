#!/usr/bin/env python
"""
    Calculates the reverse complement of a sequence.

    Usage:
    reverse.py <sequence file> <output file> <sequence file format> <output file format> <alphabet: dna or rna>
"""
import sys

from Bio import SeqIO
from Bio.Alphabet import DNAAlphabet, RNAAlphabet
from Bio.SeqRecord import SeqRecord


def main(sequence_file, ofile, iformat, oformat, alphabet="dna"):

    handle = open(sequence_file, "rU")
    recs = _reverse_complement(SeqIO.parse(handle, iformat), alphabet)
    SeqIO.write(recs, ofile, oformat)
    handle.close()


def reverse_complement(seq_iter, alphabet):
    """
    calculate the reverse complement for given sequence records
    """
    for record in seq_iter:
        if alphabet == "rna":
            record.seq.alphabet = RNAAlphabet()
        else:
            record.seq.alphabet = DNAAlphabet()
        yield SeqRecord(
            record.seq.reverse_complement(),
            description=record.description.strip(),
            id=record.id,
            name=record.name,
        )


if __name__ == "__main__":
    main(*sys.argv[1:])
