#!/usr/bin/env python

"""
    Input: molecular input file.
    Output: Molecule file with removed ions and fragments.
    Copyright 2012, Bjoern Gruening and Xavier Lucas
"""
import argparse

from openbabel import openbabel, pybel
openbabel.obErrorLog.StopLogging()


def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-iformat', default='sdf', help='input file format')
    parser.add_argument('-i', '--input', required=True, help='input file name')
    parser.add_argument('-o', '--output', required=True, help='output file name')
    parser.add_argument('-idx', default=False, action='store_true', help='should output be an indexed text table? works only for inchi/smiles, otherwise is ignored')
    return parser.parse_args()


def remove_ions(args):
    outfile = open(args.output, 'w')

    for index, mol in enumerate(pybel.readfile(args.iformat, args.input)):
        if mol.OBMol.NumHvyAtoms() > 5:
            mol.OBMol.StripSalts(0)
            if 'inchi' in mol.data:
                del mol.data['inchi']  # remove inchi cache so modified mol is saved

        mol = mol.write(args.iformat) if mol.OBMol.NumHvyAtoms() > 5 else '\n'

        if args.idx and args.iformat in ['inchi', 'smi']:
            outfile.write(f'{index}\t{mol}')
        elif mol != '\n':
            outfile.write(f'{mol}')


def __main__():
    """
        Remove any counterion and delete any fragment but the largest one for each molecule.
    """
    args = parse_command_line()
    remove_ions(args)


if __name__ == "__main__":
    __main__()
