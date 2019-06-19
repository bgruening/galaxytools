#!/usr/bin/env python
"""
Basic SuCOS scoring. Allows a set of molecules from a SD file to be ovelayed to a reference molecule,
with the resulting scores being written as properties in the output SD file.

SuCOS is the work of Susan Leung.
Bitbucket: https://bitbucket.org/Susanhleung/sucos/
Publication: https://doi.org/10.26434/chemrxiv.8100203.v1
"""

from __future__ import print_function
import argparse, os, sys, gzip
import numpy as np
from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps
import utils


### start function definitions #########################################

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
    if not result:
        utils.log("Filtered out feature type", f.GetFamily())
    return result

def getRawFeatures(mol):

    rawFeats = fdef.GetFeaturesForMol(mol)
    # filter that list down to only include the ones we're interested in
    filtered = list(filter(filterFeature, rawFeats))
    return filtered

def get_FeatureMapScore(small_feats, large_feats, tani=False):
    """
    Generate the feature map score.

    :param small_feats:
    :param large_feats:
    :param tani:
    :return:
    """

    featLists = []
    for rawFeats in [small_feats, large_feats]:
        # filter that list down to only include the ones we're intereted in
        featLists.append(rawFeats)
    fms = [FeatMaps.FeatMap(feats=x, weights=[1] * len(x), params=fmParams) for x in featLists]

    try:
        if tani:
            c = fms[0].ScoreFeats(featLists[1])
            A = fms[0].GetNumFeatures()
            B = len(featLists[1])
            if B != fms[1].GetNumFeatures():
                utils.log("Why isn't B equal to number of features...?!")
            tani_score = float(c) / (A+B-c)
            return tani_score
        else:
            fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
            return fm_score
    except ZeroDivisionError:
        utils.log("ZeroDivisionError")
        return 0

    if tani:
        tani_score = float(c) / (A+B-c)
        return tani_score
    else:
        fm_score = fms[0].ScoreFeats(featLists[1]) / min(fms[0].GetNumFeatures(), len(featLists[1]))
        return fm_score


def get_SucosScore(ref_mol, query_mol, tani=False, ref_features=None, query_features=None):
    """
    This is the key function that calculates the SuCOS scores and is expected to be called from other modules.
    To improve performance you can pre-calculate the features and pass them in as optional parameters to avoid having
    to recalculate them. Use the getRawFeatures function to pre-calculate the features.

    :param ref_mol: The reference molecule to compare to
    :param query_mol: The molecule to align to the reference
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

    fm_score = get_FeatureMapScore(ref_features, query_features, tani)
    fm_score = np.clip(fm_score, 0, 1)

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

def process(refmol_filename, inputs_filename, outputs_filename, refmol_index=None, refmol_format=None, tani=False):
    ref_mol = utils.read_single_molecule(refmol_filename, index=refmol_index, format=refmol_format)
    #utils.log("Reference mol has", ref_mol.GetNumHeavyAtoms(), "heavy atoms")
    ref_features = getRawFeatures(ref_mol)

    input_file = utils.open_file_for_reading(inputs_filename)
    suppl = Chem.ForwardSDMolSupplier(input_file)
    output_file = utils.open_file_for_writing(outputs_filename)
    writer = Chem.SDWriter(output_file)

    count = 0
    total = 0
    errors = 0
    for mol in suppl:
        count +=1
        if mol is None:
            continue
        #utils.log("Mol has", str(mol.GetNumHeavyAtoms()), "heavy atoms")
        try:
            sucos_score, fm_score, val3 = get_SucosScore(ref_mol, mol, tani=tani, ref_features=ref_features)
            mol.SetDoubleProp("SuCOS_Score", sucos_score)
            mol.SetDoubleProp("SuCOS_FeatureMap_Score", fm_score)
            if is_tanimoto:
                mol.SetDoubleProp("SuCOS_Tanimoto_Score", val3)
            else:
                mol.SetDoubleProp("SuCOS_Protrude_Score", val3)
            utils.log("Scores:", sucos_score, fm_score, val3)
            writer.write(mol)
            total +=1
        except ValueError as e:
            errors +=1
            utils.log("Molecule", count, "failed to score:", e.message)

    input_file.close()
    writer.flush()
    writer.close()
    output_file.close()

    utils.log("Completed.", total, "processed, ", count, "succeeded, ", errors, "errors")

### start main execution #########################################

def main():

    parser = argparse.ArgumentParser(description='SuCOS with RDKit')
    parser.add_argument('-i', '--input', help='Input file in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('-r', '--refmol', help='Molecule to compare against in Molfile (.mol) or SDF (.sdf) format')
    parser.add_argument('--refmol-format', help="Format for the reference molecule (mol or sdf). " +
                                                "Only needed if files don't have the expected extensions")
    parser.add_argument('--refmolidx', help='Reference molecule index in SD file if not the first', type=int, default=1)
    parser.add_argument('-o', '--output', help='Output file in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('--tanimoto', action='store_true', help='Include Tanimoto distance in score')

    args = parser.parse_args()
    utils.log("SuCOS Args: ", args)

    process(args.refmol, args.input, args.output, refmol_index=args.refmolidx, refmol_format=args.refmol_format, tani=args.tanimoto)


if __name__ == "__main__":
    main()
