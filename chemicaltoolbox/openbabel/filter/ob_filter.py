#!/usr/bin/env python
"""
    Input: set of molecules with pre-calculated physico-chemical properties
    Output: set of molecules that pass all the filters
    Copyright 2012, Bjoern Gruening and Xavier Lucas

    TODO: AND/OR conditions?
"""
import sys, os
import argparse
import cheminfolib
import json
import pybel
import shlex, subprocess

cheminfolib.pybel_stop_logging()

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', help='Input file name')
    parser.add_argument('-iformat', help='Input file format')
    parser.add_argument('-oformat', 
        default='smi',
        help='Output file format')
    parser.add_argument('-o', '--output', 
        help='Output file name',
        required=True)
    parser.add_argument('--filters', 
        help="Specify the filters to apply",
        required=True,
        )
    return parser.parse_args()

def filter_precalculated_compounds(args, filters):
    outfile = pybel.Outputfile(args.oformat, args.output, overwrite=True)
    for mol in pybel.readfile('sdf', args.input):
        for key, elem in filters.items():
            # map the short description to the larger metadata names stored in the sdf file
            property = cheminfolib.ColumnNames[key]
            min = elem[0]
            max = elem[1]
            if float(mol.data[property]) >= float(min) and float(mol.data[property]) <= float(max):
                pass
            else:
                # leave the filter loop, because one filter constrained are not satisfied
                break
        else:
            # if the filter loop terminates in a normal way (no break) all filter rules are satisfied, so save the compound
            outfile.write(mol)
    outfile.close()

def filter_new_compounds(args, filters):

    if args.iformat == args.oformat:
        # use the -ocopy option from openbabel to speed up the filtering, additionally no conversion is carried out
        # http://openbabel.org/docs/dev/FileFormats/Copy_raw_text.html#copy-raw-text
        cmd = 'obabel -i%s %s -ocopy -O %s --filter' % (args.iformat, args.input, args.output)
    else:
        cmd = 'obabel -i%s %s -o%s -O %s --filter' % (args.iformat, args.input, args.oformat, args.output)
    filter_cmd = ''
    # OBDescriptor stores a mapping from our desc shortcut to the OB name [0] and a long description [1]
    for key, elem in filters.items():
        ob_descriptor_name = cheminfolib.OBDescriptor[key][0]
        min = elem[0]
        max = elem[1]
        filter_cmd += ' %s>=%s %s<=%s ' % (ob_descriptor_name, min, ob_descriptor_name, max)

    args = shlex.split('%s "%s"' % (cmd, filter_cmd))
    #print '%s "%s"' % (cmd, filter_cmd)
    # calling openbabel with subprocess and pipe potential errors occuring in openbabel to stdout
    child = subprocess.Popen(args,
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = child.communicate()
    return_code = child.returncode

    if return_code:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        sys.stderr.write("Return error code %i from command:\n" % return_code)
        sys.stderr.write("%s\n" % cmd)
    else:
        sys.stdout.write(stdout)
        sys.stdout.write(stderr)


def __main__():
    """
        Select compounds with certain properties from a small library
    """
    args = parse_command_line()
    # Its a small trick to get the parameters in an easy way from the xml file.
    # To keep it readable in the xml file, many white-spaces are included in that string it needs to be removed.
    # Also the last loop creates a ',{' that is not an valid jason expression.
    filters = json.loads((args.filters).replace(' ', '').replace(',}', '}'))
    if args.iformat == 'sdf':
        # Check if the sdf file contains all of the required metadata to invoke the precalculation filtering
        mol = pybel.readfile('sdf', args.input).next()
        for key, elem in filters.items():
            property = cheminfolib.ColumnNames[key]
            if not property in mol.data:
                break
        else:
            # if the for loop finishes in a normal way, we should habe all properties at least in the first molecule
            # assume it is the same for all other molecules and start the precalculated filtering
            filter_precalculated_compounds(args, filters)
            return True
    filter_new_compounds(args, filters)


if __name__ == "__main__" :
    __main__()
