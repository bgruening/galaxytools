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
    return metabolic_tree.to_list()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True, help="Path to the input file.")
    parser.add_argument("-o", "--outfile", required=True, help="Path to the output file.")
    parser.add_argument("--iformat", required=True, help="Specify the input file format.")
    parser.add_argument("--phase1", required=True, help="Number of phase1 cycles.")
    parser.add_argument("--phase2", required=True, help="Number of phase2 cycles.")
    parser.add_argument("--detailed", dest="detailed",
        action="store_true", help="Returns more detailed output")
    args = parser.parse_args()

    mols = mol_supplier(args.infile, args.iformat)
    if args.detailed:
        outp = np.zeros((0,6))
    else:
        outp = np.zeros((0,3))
    for n in range(len(mols)):
        metabs = predict_metabolites(mols[n], args.phase1, args.phase2)
        for entry in range(len(metabs)):
            smiles = Chem.MolToSmiles(metabs[entry]['SyGMa_metabolite'])
            if args.detailed:
                out = np.column_stack((
                    'SYGMA{}MOL{}'.format(n, entry), # SMILES label
                    smiles, # SMILES
                    np.round(np.array(metabs[entry]['SyGMa_score'], dtype=float),
                        decimals=5), # score rounded to 5 dp
                    Chem.rdMolDescriptors.CalcMolFormula(Chem.MolFromSmiles(smiles)), # Molecular formula
                    len(metabs[entry]["SyGMa_pathway"].split("\n")), # SyGMa_n Sygma pathway length
                    metabs[entry]["SyGMa_pathway"].replace("\n", "") # SyGMa pathway
                ))
            else:
                out = np.column_stack((
                    'SYGMA{}MOL{}'.format(n, entry), # SMILES label
                    smiles, # SMILES
                    np.round(np.array(metabs[entry]['SyGMa_score'], dtype=float),
                        decimals=5) # score rounded to 5 dp
                ))
            outp = np.vstack((outp, out))
    np.savetxt(args.outfile, outp, fmt="%s", delimiter="\t")


if __name__ == "__main__":
    main()
