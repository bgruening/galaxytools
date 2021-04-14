#!/usr/bin/env python

"""
    Converts a SD-file to a GSPAN file.
"""

import os
import sys
import argparse
import openbabel
import pybel


def main(args):

    for infile in args.infile:
        file_extension = args.format or os.path.splitext(infile)[-1].lstrip(".")

        if not args.format and file_extension not in ["smi", "sdf", "inchi", "mol"]:
            sys.exit(
                "Could not guess the format from the file extension please specify with the --format option."
            )

        molecules = pybel.readfile(file_extension, infile)
        for mol in molecules:
            args.outfile.write("t # id %s\n" % mol.title.strip())
            for atom in openbabel.OBMolAtomIter(mol.OBMol):
                label = atom.GetAtomicNum()
                vertex_index = atom.GetIdx()
                args.outfile.write("v %s %s\n" % (vertex_index, label))

            for bond in openbabel.OBMolBondIter(mol.OBMol):
                src_index = bond.GetBeginAtomIdx()
                dest_index = bond.GetEndAtomIdx()
                assert src_index > 0
                assert dest_index > 0
                if bond.IsAromatic():
                    label = "a"
                elif bond.IsSingle():
                    label = "s"
                elif bond.IsDouble():
                    label = "d"
                elif bond.IsTriple():
                    label = "t"
                atom1 = bond.GetBeginAtom()
                atom2 = bond.GetEndAtom()
                args.outfile.write("e %s %s %s\n" % (src_index, dest_index, label))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-i", "--infile", nargs="*", help="Specify one or more input files"
    )
    parser.add_argument("-f", "--format", help="Format of the input file.")
    parser.add_argument(
        "--outfile",
        type=argparse.FileType("w"),
        default=sys.stdout,
        help="Specify one output file",
    )
    args = parser.parse_args()
    main(args)
