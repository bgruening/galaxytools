#!/usr/bin/python
# Script: gff_to_gb.py
# Author: Daniel Desiro'
"""
Description:
    Convert a GFF and corresponding fasta file into GenBank format.

Usage:
    gff_to_gb.py -g <gff_file> -f <fasta_file> -b <genbank_file>

Source:
    https://github.com/desiro/gffDB/blob/master/gff_to_gb.py

"""
import sys
import argparse
from Bio import SeqIO
from Bio import Seq
from BCBio import GFF
from Bio.Alphabet import generic_dna
from utils import check_gff, fix_ncbi_id


## main function
def main(gffFile, genBankFile, fastaFile):
    fastaParsed = SeqIO.to_dict(SeqIO.parse(fastaFile, "fasta", generic_dna))
    gffParsed = GFF.parse(gffFile, fastaParsed)
    SeqIO.write(check_gff(fix_ncbi_id(gffParsed)), genBankFile, "genbank")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'gff_to_gb.py', description = 'Convert a GFF and corresponding fasta file into GenBank format.', prefix_chars='-+', epilog="")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--gff', '-g', dest='gffFile', required=True, help='GFF file')
    parser.add_argument('--fasta', '-f', dest='fastaFile', required=True, help='fasta file')
    parser.add_argument('--genbank', '-b', dest='genBankFile', required=True, help='GenBank file name')
    
    options = parser.parse_args()
    main(options.gffFile, options.genBankFile, options.fastaFile)

