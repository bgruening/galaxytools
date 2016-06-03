#!/usr/bin/env python
"""
    Input:  Molecule file
    Output: Molecule file with hydrogen atoms added at the target pH.
"""
import sys, os
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel

def parse_command_line(argv):
    parser = argparse.ArgumentParser()
    parser.add_argument('--iformat', type=str, default='sdf' , help='input file format')
    parser.add_argument('-i', '--input', type=str, required=True, help='input file name')
    parser.add_argument('-o', '--output', type=str, required=True, help='output file name')
    parser.add_argument('--polar', action="store_true", default=False, help='Add hydrogen atoms only to polar atoms')
    parser.add_argument('--pH', type=float, default="7.4", help='Specify target pH value')
    return parser.parse_args()

def addh(args):
    outfile = pybel.Outputfile(args.iformat, args.output, overwrite=True)
    for mol in pybel.readfile(args.iformat, args.input):
        if mol.OBMol.NumHvyAtoms() > 5:
            mol.removeh()
            mol.OBMol.AddHydrogens(args.polar, True, args.pH)
            outfile.write(mol)
    outfile.close()

def __main__():
    """
        Add hydrogen atoms at a certain pH value
    """
    args = parse_command_line(sys.argv)
    addh(args)

if __name__ == "__main__" :
    __main__()
