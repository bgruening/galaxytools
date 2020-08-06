#!/usr/bin/env python
"""
Basic SuCOS scoring. Allows a set of molecules from a SD file to be overlayed to a reference molecule,
with the resulting scores being written as properties in the output SD file.
SuCOS is the work of Susan Leung.
GitHub: https://github.com/susanhleung/SuCOS
Publication: https://doi.org/10.26434/chemrxiv.8100203.v1
"""

from __future__ import print_function
import argparse, os, sys
import numpy as np
from rdkit import rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
from pipelines_utils import parameter_utils, utils
from pipelines_utils_rdkit import rdkit_utils

### start function definitions #########################################

field_SuCOS_Score = "SuCOS_Score"
field_SuCOS_FMScore = "SuCOS_FeatureMap_Score"
field_SuCOS_TaniScore = "SuCOS_Tanimoto_Score"
field_SuCOS_ProtrudeScore = "SuCOS_Protrude_Score"

# Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

fmParams = {}
for k in fdef.GetFeatureFamilies():
    fparams = FeatMaps.FeatMapParams()
    fmParams[k] = fparams

keep = ('Donor', 'Acceptor', 'NegIonizable', 'PosIonizable', 'ZnBinder',
        'Aromatic', 'Hydrophobe', 'LumpedHydrophobe')

def filterFeature(f):
    result = f.GetFamily() in keep
    # TODO - nothing ever seems to be filtered. Is this expected?
    #if not result:
    #    args.logfile.write("Filtered out feature type %s\n\n" % f.GetFamily())
    return result

def getRawFeatures(mol):

    rawFeats = fdef.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're interested in
    filtered = list(filter(filterFeature, rawFeats))
    return filtered

def get_FeatureMapScore(args, small_feats, large_feats, tani=False, score_mode=FeatMaps.FeatMapScoreMode.All):
    """
    Generate the feature map score.
    :param small_feats:
    :param large_feats:
    :param tani:
    :return:
    """

    featLists = []
    for rawFeats in [small_feats, large_feats]:
        # filter that list down to only include the ones we're interested in
        featLists.append(rawFeats)
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]
    # set the score mode
    fms[0].scoreMode = score_mode

    try:
        if tani:
            c = fms[0].ScoreFeats(featLists[1])
            A = fms[0].GetNumFeatures()
            B = len(featLists[1])
            if B != fms[1].GetNumFeatures():
                args.logfile.write("Why isn't B equal to number of features...?!\n\n")
            tani_score = float(c) / (A+B-c)
            return tani_score
        else:
            fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
            return fm_score
    except ZeroDivisionError:
        args.logfile.write("ZeroDivisionError\n\n")
        return 0

    if tani:
        tani_score = float(c) / (A+B-c)
        return tani_score
    else:
        fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
        return fm_score


def get_SucosScore(args, ref_mol, query_mol, tani=False, ref_features=None, query_features=None, score_mode=FeatMaps.FeatMapScoreMode.All):
    """
    This is the key function that calculates the SuCOS scores and is expected to be called from other modules.
    To improve performance you can pre-calculate the features and pass them in as optional parameters to avoid having
    to recalculate them. Use the getRawFeatures function to pre-calculate the features.
    :param ref_mol: The reference molecule to compare to
    :param query_mol: The molecule to compare to the reference
    :param tani: Whether to calculate Tanimoto distances
    :param ref_features: An optional feature map for the reference molecule, avoiding the need to re-calculate it.
    :param query_features: An optional feature map for the query molecule, avoiding the need to re-calculate it.
    :return: A tuple of 3 values. 1 the sucos score, 2 the feature map score,
        3 the Tanimoto distance or 1 minus the protrude distance
    """

    if not ref_features:
        ref_features = getRawFeatures(ref_mol)
    if not query_features:
        query_features = getRawFeatures(query_mol)

    fm_score = get_FeatureMapScore(args, ref_features, query_features, tani, score_mode)
    fm_score = np.clip(fm_score, 0, 1)

    try :
        if tani:
            tani_sim = 1 - float(rdShapeHelpers.ShapeTanimotoDist(ref_mol, query_mol))
            tani_sim = np.clip(tani_sim, 0, 1)
            SuCOS_score = 0.5*fm_score + 0.5*tani_sim
            return SuCOS_score, fm_score, tani_sim
        else:
            protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref_mol, query_mol, allowReordering=False)
            protrude_dist = np.clip(protrude_dist, 0, 1)
            protrude_val = 1.0 - protrude_dist
            SuCOS_score = 0.5 * fm_score + 0.5 * protrude_val
            return SuCOS_score, fm_score, protrude_val
    except:
        args.logfile.write("Failed to calculate SuCOS scores. Returning 0,0,0\n\n")
        return 0, 0, 0


def process(args, target_mol, inputs_supplr, writer, tani=False, score_mode=FeatMaps.FeatMapScoreMode.All):

    ref_features = getRawFeatures(target_mol)

    count = 0
    total = 0
    errors = 0
    for mol in inputs_supplr:
        count +=1
        if mol is None:
            errors +=1
            continue
        args.logfile.write("Mol has %d heavy atoms\n\n" % mol.GetNumHeavyAtoms())
        try:
            sucos_score, fm_score, val3 = get_SucosScore(args, target_mol, mol, tani=tani, ref_features=ref_features, score_mode=score_mode)
            mol.SetDoubleProp(field_SuCOS_Score, sucos_score)
            mol.SetDoubleProp(field_SuCOS_FMScore, fm_score)
            if tani:
                mol.SetDoubleProp(field_SuCOS_TaniScore, val3)
            else:
                mol.SetDoubleProp(field_SuCOS_ProtrudeScore, val3)
            args.logfile.write("Scores: %f %f %f\n\n" % (sucos_score, fm_score, val3))
            writer.write(mol)
            total += 1
        except ValueError as e:
            errors +=1
            args.logfile.write("Molecule %d failed to score: %s\n\n" % (count, e.message))

    args.logfile.write("Completed %d, processed %d, succeeded %d\n\n" % (count, total, errors))

    return count, total, errors

def parse_score_mode(value):
    if value == None or value == 'all':
        return FeatMaps.FeatMapScoreMode.All
    elif value == 'closest':
        return FeatMaps.FeatMapScoreMode.Closest
    elif value == 'best':
        return FeatMaps.FeatMapScoreMode.Best
    else:
        raise ValueError(value + " is not a valid scoring mode option")


### start main execution #########################################

def main():

    parser = argparse.ArgumentParser(description='SuCOS with RDKit')
    parameter_utils.add_default_io_args(parser)
    parser.add_argument('-r', '--refmol', help='Molecule to compare against in Molfile (.mol) or SDF (.sdf) format')
    parser.add_argument('-tm', '--target', help='Target molecule to compare against')
    parser.add_argument('-tf', '--target-format', help='Target molecule format')
    parser.add_argument('-ti', '--targetidx', help='Target molecule index in file if not the first', type=int, default=1)

    parser.add_argument('--tanimoto', action='store_true', help='Include Tanimoto distance in score')
    parser.add_argument('--score_mode', choices=['all', 'closest', 'best'],
                        help="choose the scoring mode for the feature map, default is 'all'.")
    parser.add_argument('-l', '--logfile', type=argparse.FileType('w+'),
                        default=sys.stdout, help="Path to the log file, default it sdtout")

    args = parser.parse_args()
    args.logfile.write("SuCOS Args: %s\n\n" % str(args))

    score_mode = parse_score_mode(args.score_mode)

    target_mol = rdkit_utils.read_single_molecule(args.target, index=args.targetidx, format=args.target_format)
    args.logfile.write("Target mol has %d heavy atoms\n\n" % target_mol.GetNumHeavyAtoms())

    source = "sucos.py"
    datasetMetaProps = {"source":source, "description": "SuCOS using RDKit " + rdBase.rdkitVersion}

    clsMappings = {}
    fieldMetaProps = []
    clsMappings[field_SuCOS_Score] = "java.lang.Float"
    clsMappings[field_SuCOS_FMScore] = "java.lang.Float"
    fieldMetaProps.append({"fieldName":field_SuCOS_Score,   "values": {"source":source, "description":"SuCOS score"}})
    fieldMetaProps.append({"fieldName":field_SuCOS_FMScore,   "values": {"source":source, "description":"SuCOS Feature Map score"}})

    if args.tanimoto:
        clsMappings[field_SuCOS_TaniScore] = "java.lang.Float"
        fieldMetaProps.append({"fieldName":field_SuCOS_TaniScore,   "values": {"source":source, "description":"SuCOS Tanimoto score"}})
    else:
        clsMappings[field_SuCOS_ProtrudeScore] = "java.lang.Float"
        fieldMetaProps.append({"fieldName":field_SuCOS_ProtrudeScore,   "values": {"source":source, "description":"SuCOS Protrude score"}})


    inputs_file,output,inputs_supplr,writer,output_base = rdkit_utils. \
        default_open_input_output(args.input, args.informat, args.output,
                                  'sucos', args.outformat,
                                  valueClassMappings=clsMappings,
                                  datasetMetaProps=datasetMetaProps,
                                  fieldMetaProps=fieldMetaProps)

    # this does the processing
    count, total, errors = process(args, target_mol, inputs_supplr, writer, tani=args.tanimoto, score_mode=score_mode)

    inputs_file.close()
    writer.flush()
    writer.close()
    output.close()

    if args.meta:
        utils.write_metrics(output_base, {'__InputCount__':count, '__OutputCount__':total, '__ErrorCount__':errors, 'RDKitSuCOS':total})


if __name__ == "__main__":
    main()