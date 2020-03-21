#!/usr/bin/env python3
import argparse
import pandas as pd
from rdkit import Chem

def sdf_to_tab(vars):
    mols = Chem.SDMolSupplier(vars.inp, sanitize=False)
    df = pd.DataFrame()  # for output

    for n in range(len(mols)):
        if mols[n]:
            d = mols[n].GetPropsAsDict()
            # filter dict for desired props
            if vars.props.strip() == '':  # none specified, return all
                d = {prop: val for (prop, val) in d.items() if not any(x in str(val) for x in ['\n', '\t'])}  # remove items containing newlines or tabs
            else:
                d = {prop: val for (prop, val) in d.items() if prop in vars.props.replace(' ', '').split(',')}  # remove items not requested via CLI
            if vars.name:
                d['SDFMoleculeName'] = mols[n].GetProp('_Name')
            if vars.smiles:
                d['SMILES'] = Chem.MolToSmiles(mols[n], isomericSmiles=False)
            d['Index'] = int(n)

            df = df.append(d, ignore_index=True)
        else:
            print("Molecule could not be read - skipped.")

    df = df.astype({'Index': int}).set_index('Index')
    df.to_csv(vars.out, sep='\t', header=vars.header)

def main():
    parser = argparse.ArgumentParser(description="Convert SDF to tabular")
    parser.add_argument('--inp', '-i', help="The input file", required=True)
    parser.add_argument('--out', '-o', help="The output file", required=True)
    parser.add_argument('--props', '-p', help="Properties to filter (leave blank for all)", required=True)
    parser.add_argument('--header', '-t', action='store_true',
                        help="Write property name as the first row.")
    parser.add_argument('--smiles', '-s', action='store_true',
                        help="Include SMILES in output.")
    parser.add_argument('--name', '-n', action='store_true',
                        help="Include molecule name in output.")
    sdf_to_tab(parser.parse_args())
    

if __name__ == "__main__":
    main()
