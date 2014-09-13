#!/usr/bin/env python
"""
    Input: Molecular input file.
    Output: Physico-chemical properties are computed and stored as metadata in the sdf output file.
    Copyright 2012, Bjoern Gruening and Xavier Lucas
"""
import sys, os
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel
import cheminfolib


def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--iformat', default='sdf' , help='input file format')
    parser.add_argument('-i', '--input', required=True, help='input file name')
    parser.add_argument('--oformat', default='sdf', choices = ['sdf', 'table'] , help='output file format')
    parser.add_argument('--header', type=bool, help='Include the header as the first line of the output table')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    return parser.parse_args()

def compute_properties(args):
    if args.oformat == 'sdf':
        outfile = pybel.Outputfile(args.oformat, args.output, overwrite=True)
    else:
        outfile = open(args.output, 'w')
        if args.header:
            mol = pybel.readfile(args.iformat, args.input).next()
            metadata = cheminfolib.get_properties_ext(mol)
            outfile.write( '%s\n' % '\t'.join( [ cheminfolib.ColumnNames[key] for key in metadata ] ) )

    for mol in pybel.readfile(args.iformat, args.input):
        if mol.OBMol.NumHvyAtoms() > 5:
            metadata = cheminfolib.get_properties_ext(mol)
            if args.oformat == 'sdf':
                [ mol.data.update( { cheminfolib.ColumnNames[key] : metadata[key] } ) for key in metadata ]
                outfile.write(mol)
            else:
                outfile.write( '%s\n' % ('\t'.join( [ str(metadata[key]) for key in metadata ] ) ) )
    outfile.close()

def __main__():
    """
        Physico-chemical properties are computed and stored as metadata in the sdf output file
    """
    args = parse_command_line(sys.argv)
    compute_properties(args)

if __name__ == "__main__" :
    __main__()
