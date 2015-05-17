#!/usr/bin/env python
"""
Input: UniProt accession number(s)
Output: uniprot protein sequence(s)
Sept 2011
Xavier Lucas, Bjoern Gruening
"""
import sys, os
import argparse
import psycopg2.extras
import subprocess

LIB_PATH = '/media/data/databases/uniprot/devided/'

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', required=False, type=str, help='input file name containing one or multiple unique molecule identifiers')
    parser.add_argument('-t', '--text', required=False, type=str, help='todo')
    parser.add_argument('-c', '--column', required=False, type=int, help='#column containing the id codes')
    parser.add_argument('-o', '--output', type=str, help='output file name')
    return parser.parse_args()

def fetch_fasta( identifiers, args ):
    for elem in identifiers:
        path = '%s%s/%s.fasta' % ( LIB_PATH, elem[:2], elem )
        cmd = 'cat %s >> %s' % ( path, args.output )
        try:
            subprocess.call( cmd, shell=True )
        except:
            print 'UniProt entry %s not found' % elem
            continue

def __main__():
    """
        Fetch protein sequences by UniProt entry
    """
    
    args = parse_command_line(sys.argv)

    identifiers = []
    if args.input and args.column:
        [ identifiers.append( line.split('\t')[args.column - 1].strip() ) for line in open(args.input, 'r') ]
    elif args.text:
        identifiers = [ a.strip() for a in args.text.split() if a.strip() ]

    fetch_fasta( identifiers, args )

if __name__ == "__main__" :
    __main__()
