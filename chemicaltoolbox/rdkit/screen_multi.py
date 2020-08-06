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

from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem.Fingerprints import FingerprintMols

from pipelines_utils_rdkit import filter, rdkit_utils, mol_utils
from pipelines_utils import parameter_utils, utils

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


def main():
    ### command line args defintions #########################################

    parser = argparse.ArgumentParser(description='RDKit screen')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('--qsmiles',
                       help='filename of query structures as smiles (incompatible with --sdf and --qjson args)')
    group.add_argument('--qsdf',
                       help='filename of query structures as sdfile (incompatible with --smiles and --qjson args)')
    group.add_argument('--qjson',
                       help='filename of query structures as MoleculeObject JSON (incompatible with --qsmiles and --qsdf args)')
    parser.add_argument('--qsmilesTitleLine', action='store_true', help='the smiles file has a title line')
    parser.add_argument('--qsmilesDelimiter', default='\t', help='delimiter for smiles file (default is tab)')
    parser.add_argument('--qsmilesColumn', type=int, default=0,
                        help='column in smiles file with the smiles (default is first column)')
    parser.add_argument('--qsmilesNameColumn', type=int, default=1,
                        help='column in smiles file with ID (default is second column)')
    parser.add_argument('--qprop',
                        help='property name in query molecules to report. If not defined (or property is not present) ' +
                             'then name property is not written. JSON format uses the UUID as default')

    parser.add_argument('--simmin', type=float, default=0.7, help='similarity lower cutoff (1.0 means identical)')
    parser.add_argument('--simmax', type=float, default=1.0, help='similarity upper cutoff (1.0 means identical)')
    parser.add_argument('-d', '--descriptor', type=str.lower, choices=list(descriptors.keys()), default='rdkit',
                        help='descriptor or fingerprint type (default rdkit)')
    parser.add_argument('-m', '--metric', type=str.lower, choices=list(metrics.keys()), default='tanimoto',
                        help='similarity metric (default tanimoto)')
    parser.add_argument('-f', '--fragment', choices=['hac', 'mw'],
                        help='Find single fragment if more than one (hac = biggest by heavy atom count, mw = biggest by mol weight )')
    parser.add_argument('--hacmin', type=int, help='Min heavy atom count')
    parser.add_argument('--hacmax', type=int, help='Max heavy atom count')
    parser.add_argument('--mwmin', type=float, help='Min mol weight')
    parser.add_argument('--mwmax', type=float, help='Max mol weight')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('--thin', action='store_true', help='Thin output mode')
    parser.add_argument('-q', '--quiet', action='store_true', help='Quiet mode')
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")

    args = parser.parse_args()
    args.logfile.write("Screen Args: %s\n\n" % str(args))

    descriptor = descriptors[args.descriptor.lower()]
    metric = metrics[args.metric.lower()]

    propName = args.qprop
    if args.qsmiles:
        queryMolsupplier = rdkit_utils.default_open_input_smiles(args.qsmiles, delimiter=args.qsmilesDelimiter,
                                                                 smilesColumn=args.qsmilesColumn,
                                                                 nameColumn=args.qsmilesNameColumn,
                                                                 titleLine=args.qsmilesTitleLine)
        queryInput = None
    elif args.qsdf:
        queryInput, queryMolsupplier = rdkit_utils.default_open_input_sdf(args.qsdf)
    elif args.qjson:
        queryInput, queryMolsupplier = rdkit_utils.default_open_input_json(args.qjson, lazy=False)
        if not propName:
            propName = "uuid"
    else:
        raise ValueError('No query structure specified')

    queryFps = {}
    args.logfile.write("Preparing query fingerprints\n\n")
    count = 0
    for q in queryMolsupplier:
        count += 1
        if q:
            queryFps[q] = descriptor(q)
        else:
            args.logfile.write("WARNING: Failed to parse Molecule %d\n\n" % count)
    if queryInput:
        queryInput.close()

    input, output, suppl, writer, output_base = rdkit_utils.default_open_input_output(args.input, args.informat,
                                                                                      args.output, 'screen_multi',
                                                                                      args.outformat)

    # OK, all looks good so we can hope that things will run OK.
    # But before we start lets write the metadata so that the results can be handled.
    # if args.meta:
    #    t = open(output_base + '_types.txt', 'w')
    #    t.write(field_Similarity + '=integer\n')
    #    t.flush()
    #    t.close()

    writer = rdkit_utils.ThickSDWriter(args.output)
    i = 0
    count = 0
    for mol in suppl:
        i += 1
        if mol is None: continue
        if args.fragment:
            mol = mol_utils.fragment(mol, args.fragment, quiet=args.quiet)
        if not filter.filter(mol, minHac=args.hacmin, maxHac=args.hacmax, minMw=args.mwmin, maxMw=args.mwmax,
                             quiet=args.quiet):
            continue
        targetFp = descriptor(mol)
        idx = 0
        hits = 0
        bestScore = 0
        bestName = None
        for queryMol in queryFps:
            idx += 1
            sim = metric(queryFps[queryMol], targetFp)
            if propName:
                name = str(queryMol.GetProp(propName))
            else:
                name = None
            if sim >= args.simmin and sim <= args.simmax:
                hits += 1
                if not args.quiet:
                    args.logfile.write("%d %d %f\n\n" % (i, idx, sim))
                if sim > bestScore:
                    bestScore = sim
                    bestIdx = idx
                    if name:
                        bestName = name
                if name:
                    mol.SetDoubleProp(field_Similarity + "_" + name, sim)
                else:
                    mol.SetDoubleProp(field_Similarity + "_" + str(idx) + "_Score", sim)

        if hits > 0:
            count += 1
            mol.SetDoubleProp(field_Similarity + "_BestScore", bestScore)
            if bestName:
                mol.SetProp(field_Similarity + "_BestName", bestName)
            else:
                mol.SetIntProp(field_Similarity + "_BestIndex", bestIdx)
            mol.SetIntProp(field_Similarity + "_Count", hits)
            writer.write(mol)

    args.logfile.write("Found %d similar molecules\n\n" % count)

    writer.flush()
    writer.close()
    input.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__': i, '__OutputCount__': count, 'RDKitScreen': count})

    return count


if __name__ == "__main__":
    main()