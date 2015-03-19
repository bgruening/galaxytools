#!/usr/bin/env python

import re
import sys
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
from Bio.Alphabet import generic_dna

choices = ['embl', 'fasta', 'fastq-sanger', 'fastq', 'fastq-solexa', 'fastq-illumina', 'genbank', 'gb']

def find_pattern(seqs, pattern, outfile_path):
    """
    Finds all occurrences of a pattern in the a given sequence.
    Outputs sequence ID, start and end postion of the pattern.
    """
    pattern = pattern.upper()
    rev_compl = Seq(pattern, generic_dna).complement()
    search_func = simple_pattern_search
    if set(pattern).difference(set('ATCG')):
        search_func = complex_pattern_search

    with open(outfile_path, 'w+') as outfile:
        for seq in seqs:
            search_func(seq, pattern, outfile)
            search_func(seq, rev_compl, outfile, '-')


def simple_pattern_search(sequence, pattern, outfile, strand='+'):
    """
    Simple regular expression search. This is way faster than the complex search.
    """
    bed_template = '%s\t%s\t%s\t%s\t%s\t%s\n'
    for match in re.finditer( str(pattern), str(sequence.seq) ):
        outfile.write(bed_template % (sequence.id,  match.start(), match.end(), sequence.name, '', strand))


def complex_pattern_search(sequence, pattern, outfile, strand='+'):
    """
    Searching for pattern with biopyhon's nt_search().
    This allows for ambiguous values, like N = A or T or C or G, R = A or G ...
    """
    l = len(pattern)
    matches = nt_search(str(sequence.seq), pattern)
    bed_template = '%s\t%s\t%s\t%s\t%s\t%s\n'
    for match in matches[1:]:
        outfile.write(bed_template % (sequence.id, match, match+l, sequence.name, '', strand) )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input' , required=True)
    parser.add_argument('-o', '--output' , required=True)
    parser.add_argument('-p', '--pattern' , required=True)
    parser.add_argument('-f', '--format', default="fasta", choices=choices)
    args = parser.parse_args()

    with open(args.input) as handle:
        find_pattern( SeqIO.parse(handle, args.format), args.pattern, args.output )

