#!/usr/bin/env python3
import argparse
import pandas as pd
from rdkit import Chem

from rdkit_util import get_supplier

def get_formal_charge(mols):
    return [Chem.GetFormalCharge(mol) for mol in mols]

def main():
    parser = argparse.ArgumentParser(description="Convert SDF to tabular")
    parser.add_argument('--infile', '-i', help="The input file", required=True)
    parser.add_argument('--outfile', '-o', help="The output file", required=True)
    parser.add_argument('--iformat', '-f', choices=['sdf', 'smiles', 'inchi', 'mol2', 'pdb'], help="Format for the input file.", required=True)
    # parser.add_argument('--smiles', '-s', action='store_true',
    #                     help="Include SMILES in output.")
    # parser.add_argument('--name', '-n', action='store_true',
    #                     help="Include molecule name in output.")

    args = parser.parse_args()

    if args.iformat in ['inchi', 'smiles']:
        mols = get_supplier(args.infile, format='inchi')
    elif args.iformat == 'sdf':
        mols = Chem.SDMolSupplier(args.infile)
    elif args.iformat == 'pdb':
        mols = [Chem.MolFromPDBFile(args.infile)]
    elif args.iformat == 'mol2':
        mols = [Chem.MolFromMol2File(args.infile)]

    charges = [str(c) for c in get_formal_charge(mols)]

    with open(args.outfile, 'w') as f:
        f.write('\n'.join(charges))
    

if __name__ == "__main__":
    main()
