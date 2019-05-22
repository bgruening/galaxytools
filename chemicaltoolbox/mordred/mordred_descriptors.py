import argparse
import numpy as np
import pandas as pd
from mordred import Calculator, descriptors
from mordred.error import Missing, Error
from rdkit import Chem
from rdkit.Chem.rdmolfiles import SDMolSupplier, SmilesMolSupplier


def convert_errors_to_nan(el):
    """
    Remove elements from the Mordred dataframe which are not
    in float or int format 
    """
    if type(el) == bool:
        return int(el)
    if type(el) not in [float, int, np.float64]:
        return None
    return el


def mol_supplier(filename, ext):
    """
    Based on the file extension, use the appropriate RDKit function to
    load a chemical data file (SMILES or SDF) containing multiple molecules
    and return a list of RDKit Mol objects
    """
    if ext == 'sdf':
        return [n for n in SDMolSupplier(filename)]
    with open(filename) as f: 
        mols = f.read().split('\n') 
    if ext == 'smi':
        return [Chem.MolFromSmiles(mol, sanitize=True) for mol in mols]
    if ext == 'inchi':
        return [Chem.inchi.MolFromInchi(mol, sanitize=True) for mol in mols]


def mordred_descriptors(mols, output, header, use_3d):
    """
    Calculate Mordred descriptors and save as tabular
    """
    calc = Calculator(descriptors, ignore_3D=(not use_3d))
    invalid_mols = np.where(np.array(mols) == None)[0]  # indices of invalid SMILES/SDMols
    mols = [Chem.MolFromSmiles('') if n is None else n for n in mols]  # replace invalid mols with placeholder
    df = calc.pandas(mols, quiet=True)  # calculate descriptors
    for mol in invalid_mols:  # remove placeholders
        df.iloc[mol] = np.nan
    df = df.applymap(convert_errors_to_nan)  # remove descriptors which errored
    df.to_csv(output, na_rep='', sep='\t', index=False, header=header)  # write output


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Path to the input file.')
    parser.add_argument("--iformat", help="Specify the input file format.")

    parser.add_argument('-o', '--outfile', type=argparse.FileType('w+'), 
                    help="Path to the result file")

    parser.add_argument("--3d", dest="use_3d", action="store_true",
                    default=False,
                    help="Use 3d descriptors - only with SDF input.")

    parser.add_argument("--header", dest="header", action="store_true",
                    default=False,
                    help="Write header line.")
    args = parser.parse_args()

    mols = mol_supplier(args.infile, args.iformat)
    mordred_descriptors(mols, args.outfile, args.header, args.use_3d)