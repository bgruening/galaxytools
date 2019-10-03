#!/usr/bin/env python3

import argparse
import csv
import sygma
import numpy as np
from rdkit import Chem
from rdkit.Chem.rdmolfiles import SDMolSupplier, SmilesMolSupplier

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
    if ext == 'smi' or ext == 'inchi':
        return [Chem.MolFromSmiles(mol, sanitize=True) for mol in mols if mol != '']

def predict_metabolites(parent, phase1_cycles, phase2_cycles):
    """
    Prediction of metabolites derived from a parent molecule
    """
    scenario = sygma.Scenario([
        [sygma.ruleset['phase1'], int(phase1_cycles)],
        [sygma.ruleset['phase2'], int(phase2_cycles)]])
    metabolic_tree = scenario.run(parent)
    metabolic_tree.calc_scores()
    return metabolic_tree.to_smiles()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Path to the input file.')
    parser.add_argument('-o', '--outfile', required=True, help='Path to the output file.')
    parser.add_argument("--iformat", help="Specify the input file format.")
    parser.add_argument("--phase1", help="Number of phase1 cycles.")
    parser.add_argument("--phase2", help="Number of phase2 cycles.")
    args = parser.parse_args()

    mols = mol_supplier(args.infile, args.iformat)
    outp = np.zeros((0,3))
    for n in range(len(mols)):
        metabs = np.array(predict_metabolites(mols[n], args.phase1, args.phase2))
        metabs = np.column_stack((
            metabs[:,0],  # SMILES
            ['SYGMA{}MOL{}'.format(n, m) for m in range(metabs.shape[0])],  # SMILES label
            np.round(np.array(metabs[:,1], dtype=float), decimals=5)  # score rounded to 5 dp
        ))
        outp = np.vstack((outp, metabs))
    np.savetxt(args.outfile, outp, fmt="%s")


if __name__ == "__main__":
    main()
