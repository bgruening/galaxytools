#!/usr/bin/env python

from rdkit.Chem import Descriptors
from rdkit import Chem
import sys, os, re
import argparse
import inspect

def get_supplier( infile, format = 'smiles' ):
    """
    Returns a generator over a SMILES or InChI file. Every element is of RDKit 
    molecule and has its original string as _Name property.
    """
    with open(infile) as handle:
        for line in handle:
            line = line.strip()
            if format == 'smiles':
                mol = Chem.MolFromSmiles( line, sanitize=True )
            elif format == 'inchi':
                mol = Chem.inchi.MolFromInchi( line, sanitize=True, removeHs=True, logLevel=None, treatWarningAsError=False )
            if mol is None:
                yield False
            else:
                mol.SetProp( '_Name', line.split('\t')[0] )
                yield mol

def get_rdkit_descriptor_functions():
    """
    Returns all descriptor functions under the Chem.Descriptors Module as tuple of (name, function)
    """
    ret = [ (name, f) for name, f in inspect.getmembers( Descriptors ) if inspect.isfunction( f ) and not name.startswith( '_' ) ]
    ret.sort()
    return ret


def descriptors( mol, functions ):
    """
    Calculates the descriptors of a given molecule.
    """
    for name, function in functions:
        yield (name, function( mol ))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Path to the input file.')
    parser.add_argument("--iformat", help="Specify the input file format.")

    parser.add_argument('-o', '--outfile', type=argparse.FileType('w+'), 
        default=sys.stdout, help="path to the result file, default it sdtout")

    parser.add_argument("--header", dest="header", action="store_true",
                    default=False,
                    help="Write header line.")

    args = parser.parse_args()

    if args.iformat == 'sdf':
        supplier = Chem.SDMolSupplier( args.infile )
    elif args.iformat =='smi':
        supplier = get_supplier( args.infile, format = 'smiles' )
    elif args.iformat == 'inchi':
        supplier = get_supplier( args.infile, format = 'inchi' )

    functions = get_rdkit_descriptor_functions()

    if args.header:
        args.outfile.write( '%s\n' % '\t'.join( [name for name, f in functions] ) )

    for mol in supplier:
        if not mol:
            continue
        descs = descriptors( mol, functions )
        molecule_id = mol.GetProp("_Name")
        args.outfile.write( "%s\n" % '\t'.join( [molecule_id]+ [str(res) for name, res in descs] ) )

