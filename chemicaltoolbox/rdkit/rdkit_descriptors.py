#!/usr/bin/env python

import argparse
import inspect
import sys

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit_util import get_supplier


def get_rdkit_descriptor_functions():
    """
    Returns all descriptor functions under the Chem.Descriptors Module as tuple of (name, function)
    """
    ret = [(name, f) for name, f in inspect.getmembers(Descriptors) if inspect.isfunction(f) and not name.startswith('_')]
    # some which are not in the official Descriptors module we need to add manually
    ret.extend([('FormalCharge', Chem.GetFormalCharge), ('SSSR', Chem.GetSSSR)])
    ret.sort()
    return ret


def descriptors(mol, functions):
    """
    Calculates the descriptors of a given molecule.
    """
    for name, function in functions:
        yield (name, function(mol))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Path to the input file.')
    parser.add_argument("--iformat", help="Specify the input file format.")

    parser.add_argument('-o', '--outfile', type=argparse.FileType('w+'),
                        default=sys.stdout,
                        help="path to the result file, default is stdout")

    parser.add_argument('-s', '--select', default=None,
                        help="select a subset of comma-separated descriptors to use")

    parser.add_argument("--header", dest="header", action="store_true",
                        default=False,
                        help="Write header line.")

    args = parser.parse_args()

    if args.iformat == 'sdf':
        supplier = Chem.SDMolSupplier(args.infile)
    elif args.iformat == 'smi':
        supplier = get_supplier(args.infile, format='smiles')
    elif args.iformat == 'inchi':
        supplier = get_supplier(args.infile, format='inchi')
    elif args.iformat == 'pdb':
        supplier = [Chem.MolFromPDBFile(args.infile)]
    elif args.iformat == 'mol2':
        supplier = [Chem.MolFromMol2File(args.infile)]

    functions = get_rdkit_descriptor_functions()
    if args.select and args.select != 'None':
        selected = args.select.split(',')
        functions = [(name, f) for name, f in functions if name in selected]

    if args.header:
        args.outfile.write('%s\n' % '\t'.join(['MoleculeID'] + [name for name, f in functions]))

    for mol in supplier:
        if not mol:
            continue
        descs = descriptors(mol, functions)
        try:
            molecule_id = mol.GetProp("_Name")
        except KeyError:
            molecule_id = Chem.MolToSmiles(mol)
        args.outfile.write("%s\n" % '\t'.join([molecule_id] + [str(round(res, 6)) for name, res in descs]))
