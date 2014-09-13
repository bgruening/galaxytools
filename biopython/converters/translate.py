#!/usr/bin/env python
"""
    Translates a nucleotide sequences into amino acid sequences.

    Usage:
    translate.py <sequence file> <output file> <sequence file format> <output file format> <codon table number> <is complete cds? true or false> <quiet at stop codons: true or false>
"""
import sys, os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def main(sequence_file, ofile, iformat, oformat, codon_table, is_complete_cds = False, stop = False):

    if is_complete_cds == 'true':
        is_complete_cds = True
    else:
        is_complete_cds = False

    if stop == "true":
        stop = True
    else:
        stop = False

    handle = open(sequence_file, "rU")
    recs = translate(SeqIO.parse(handle, iformat), codon_table, is_complete_cds, stop)
    SeqIO.write(recs, ofile, oformat)
    handle.close()


def translate(seq_iter, codon_table, is_complete_cds, stop):
    """
        calculate the reverse complement for given sequence records
    """
    for record in seq_iter:
        yield SeqRecord( record.seq.translate(table = codon_table, to_stop = stop, cds = is_complete_cds),
            description = record.description.strip(), 
            id = record.id, 
            name = record.name )

if __name__ == "__main__" :
    main(*sys.argv[1:])
