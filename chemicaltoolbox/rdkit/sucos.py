#!/usr/bin/env python

# Copyright 2018 Informatics Matters Ltd.
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

import argparse, os, sys

from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps

from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

### start field name defintions #########################################

field_SuCOS_Score = "SuCOS_Score"

### start function definitions #########################################

#################################################
#### Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))
#    keep = ('Donor','Acceptor','NegIonizable','PosIonizable','Aromatic')

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

def get_FeatureMapScore(ref_mol, query_mol):
    featLists = []
    for m in [ref_mol, query_mol]:
        rawFeats = fdef.GetFeaturesForMol(m)
        # filter that list down to only include the ones we're interested in
        featLists.append([f for f in rawFeats if f.GetFamily() in keep])
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    fm_score = fms[0].ScoreFeats(featLists[1]) / float(min(fms[0].GetNumFeatures(), len(featLists[1])))
    return fm_score

def get_SucosScore(ref_mol, query_mol, field_name):
    fm_score = get_FeatureMapScore(ref_mol, query_mol)
    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref_mol, query_mol, allowReordering=False)
    score = 0.5 * fm_score + 0.5 * (1.0 - protrude_dist)
    query_mol.SetDoubleProp(field_name, score)
    return score

### start main execution #########################################

def main():

    parser = argparse.ArgumentParser(description='SuCOS with RDKit')
    parser.add_argument('--target', help='molecule to compare against')
    parser.add_argument('--targetidx', help="Target molecule index in SD file if not the first", type=int, default=1)
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")
    parameter_utils.add_default_io_args(parser)

    args = parser.parse_args()
    args.logfile.write("SuCOS Args: %s\n\n" % str(args))

    # TODO - handle molecules with multiple fragments

    ref_mol = rdkit_utils.read_single_molecule(args.target, index=args.targetidx)
    args.logfile.write("Reference mol has %d heavy atoms\n" % ref_mol.GetNumHeavyAtoms())

    source = "sucos.py"
    datasetMetaProps = {"source":source, "description": "SuCOS using RDKit " + rdBase.rdkitVersion}
    clsMappings = { "SuCOS_score":   "java.lang.Float" }
    fieldMetaProps = [
        {"fieldName":field_SuCOS_Score,   "values": {"source":source, "description":"SuCOS score"}}
    ]

    input,output,suppl,writer,output_base = rdkit_utils.\
        default_open_input_output(args.input, args.informat, args.output,
                                  'sucos', args.outformat,
                                  valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)

    count = 0
    total = 0
    errors = 0
    writer = rdkit_utils.ThickSDWriter(args.output)

    for mol in suppl:
        count +=1
        if mol is None:
            continue
        try:
            fm_score = get_SucosScore(ref_mol, mol, field_SuCOS_Score)
            args.logfile.write("Score: %.3f\n" % fm_score)
            writer.write(mol)
            total +=1
        except ValueError as e:
            errors +=1
            args.logfile.write("Molecule %d failed to score: %s\n" % (count, e.message))

    input.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':count, '__OutputCount__':total, '__ErrorCount__':errors, 'RDKitSuCOS':total})

if __name__ == "__main__":
    main()