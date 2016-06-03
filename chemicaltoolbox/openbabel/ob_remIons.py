#!/usr/bin/env python
"""
    Input: molecular input file.
    Output: Molecule file with removed ions and fragments.
    Copyright 2012, Bjoern Gruening and Xavier Lucas
"""
import sys, os
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-iformat', default='sdf' , help='input file format')
    parser.add_argument('-i', '--input', required=True, help='input file name')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    return parser.parse_args()

def remove_ions(args):
    outfile = pybel.Outputfile(args.iformat, args.output, overwrite=True)
    for mol in pybel.readfile(args.iformat, args.input):
        if mol.OBMol.NumHvyAtoms() > 5:
            mol.OBMol.StripSalts(0)
            # Check if new small fragments have been created and remove them
            if mol.OBMol.NumHvyAtoms() > 5:
                outfile.write(mol)
    outfile.close()

def __main__():
    """
        Remove any counterion and delete any fragment but the largest one for each molecule.
    """
    args = parse_command_line()
    remove_ions(args)

if __name__ == "__main__" :
    __main__()
