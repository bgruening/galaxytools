#!/usr/bin/env python
"""
    Input: molecular input file.
    Output: Molecule file with removed ions and fragments.
    Copyright 2013, Bjoern Gruening and Xavier Lucas
"""
import sys, os
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('--iformat', default='sdf' , help='input file format')
    parser.add_argument('-i', '--input', required=True, help='input file name')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    return parser.parse_args()

def remove_protonation( args ):
    outfile = pybel.Outputfile(args.iformat, args.output, overwrite=True)
    for mol in pybel.readfile(args.iformat, args.input):
        [atom.OBAtom.SetFormalCharge(0) for atom in mol.atoms]
        outfile.write( mol )
    outfile.close()

def __main__():
    """
        Remove any protonation state from each atom in each molecule.
    """
    args = parse_command_line()
    remove_protonation( args )

if __name__ == "__main__" :
    __main__()
