#!/usr/bin/env python

# Copyright 2017 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse, sys

from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

from pipelines_utils_rdkit import filter, rdkit_utils, mol_utils
from pipelines_utils import parameter_utils, utils

### start function defintions #########################################


def open_file(filename):
    with open(filename) as f:
        return f.readline()


### start field name defintions #########################################

field_Similarity = "Similarity"

### start main execution #########################################

descriptors = {
    # 'atompairs':   lambda m: Pairs.GetAtomPairFingerprint(m),
    'maccs': lambda m: MACCSkeys.GenMACCSKeys(m),
    'morgan2': lambda m: AllChem.GetMorganFingerprint(m, 2),
    'morgan3': lambda m: AllChem.GetMorganFingerprint(m, 3),
    'rdkit': lambda m: FingerprintMols.FingerprintMol(m),
    # 'topo':        lambda m: Torsions.GetTopologicalTorsionFingerprint(m)
}

metrics = {
    'asymmetric': DataStructs.AsymmetricSimilarity,
    'braunblanquet': DataStructs.BraunBlanquetSimilarity,
    'cosine': DataStructs.CosineSimilarity,
    'dice': DataStructs.DiceSimilarity,
    'kulczynski': DataStructs.KulczynskiSimilarity,
    'mcconnaughey': DataStructs.McConnaugheySimilarity,
    # 'onbit':DataStructs.OnBitSimilarity,
    'rogotgoldberg': DataStructs.RogotGoldbergSimilarity,
    'russel': DataStructs.RusselSimilarity,
    'sokal': DataStructs.SokalSimilarity,
    'tanimoto': DataStructs.TanimotoSimilarity
    # 'tversky': DataStructs.TverskySimilarity
}


def get_supplier(infile, format='smiles' ):
    """
    Returns a generator over a SMILES file. Every element is of RDKit
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


def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit screen')
    parser.add_argument('--smiles_mol', help='SMILES string or mol file input.')
    parser.add_argument('--isMol', action='store_true', help='SMILES string or mol file input.')
    parser.add_argument('--simmin', type=float, default=0.7, help='similarity lower cutoff (1.0 means identical)')
    parser.add_argument('--simmax', type=float, default=1.0, help='similarity upper cutoff (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit',
                        help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto',
                        help='similarity metric (default tanimoto)')
    parser.add_argument('-f', '--fragment', choices=['none', 'hac', 'mw'],
                        help='Find single fragment if more than one (hac = biggest\
                         by heavy atom count, mw = biggest by mol weight )')
    parser.add_argument('--hacmin', type=int, help='Min heavy atom count')
    parser.add_argument('--hacmax', type=int, help='Max heavy atom count')
    parser.add_argument('--mwmin', type=float, help='Min mol weight')
    parser.add_argument('--mwmax', type=float, help='Max mol weight')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-t', '--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
        default=sys.stdout, help="Path to the log file, default it sdtout")
    args = parser.parse_args()
    args.logfile.write("Screen Args: %s\n" % str(args))

    descriptor = descriptors[args.descriptor.lower()]
    metric = metrics[args.metric.lower()]

    if args.smiles_mol:
        if args.isMol:
            query_rdkitmol = Chem.MolFromMolFile(args.smiles_mol)
        else:
            query_rdkitmol = Chem.MolFromSmiles(args.smiles_mol)


    query_fp = descriptor(query_rdkitmol)


    input, output, suppl, writer, output_base = rdkit_utils.default_open_input_output(args.input, args.informat,
                                                                                        args.output, 'screen',
                                                                                        args.outformat,
                                                                                        thinOutput=args.thin)

    i = 0
    count = 0
    writer = rdkit_utils.ThickSDWriter(args.output)

    for mol in suppl:
        i += 1
        if mol is None: continue
        if args.fragment:
            mol = mol_utils.fragment(mol, args.fragment, quiet=args.quiet)
        if not filter.filter(mol, minHac=args.hacmin, maxHac=args.hacmax, minMw=args.mwmin, maxMw=args.mwmax,
                                quiet=args.quiet):
            continue
        target_fp = descriptor(mol)
        sim = metric(query_fp, target_fp)

        if sim >= args.simmin and sim <= args.simmax:
            count += 1
            if not args.quiet:
                args.logfile.write("%d %f\n" % (i, sim))
            mol.SetDoubleProp(field_Similarity, sim)
            writer.write(mol)

    args.logfile.write("Found %d similar molecules" % count)

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RDKitScreen': i})

    return count


if __name__ == "__main__":
    main()