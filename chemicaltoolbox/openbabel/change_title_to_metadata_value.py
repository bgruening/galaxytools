#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
    Change the title from a molecule file to metadata
    value of a given-id of the same molecule file.
"""

import os
import sys
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel
import random
import string


def main():
    parser = argparse.ArgumentParser(
        description="Change the title from a molecule file to metadata \
value of a given-id of the same molecule file.",
    )
    parser.add_argument('--infile', '-i', 
        required=True, help="path to the input file")
    parser.add_argument('--outfile', '-o', 
        required=True, help="path to the output file")
    parser.add_argument('--key', '-k',
        required=True, help="the metadata key from the sdf file which should inlcude the new title")
    parser.add_argument('--random', '-r',
        action="store_true", help="Add random suffix to the title.")

    args = parser.parse_args()

    output = pybel.Outputfile("sdf", args.outfile, overwrite=True)



    for mol in pybel.readfile("sdf", args.infile):
        if args.key in mol.data:
            mol.title = mol.data[args.key]
            if args.random:
                print('fooo')
                suffix = ''.join(random.choice(string.ascii_lowercase + string.digits) for _ in range(13))
                mol.title += '__%s' % suffix
        output.write( mol )




    output.close()


if __name__ == "__main__":
    main()

