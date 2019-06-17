#!/usr/bin/env python

# Original code from Tim Dudgeon<tdudgeon@informaticsmatters.com> adapted from the SuCOS service
# from the Squonk Computational Notebook:
# https://github.com/InformaticsMatters/pipelines/blob/master/src/python/pipelines/rdkit/sucos.py
#
# SuCOS is the work of Susan Leung.

from __future__ import print_function
import argparse, os, sys, gzip
from rdkit import Chem, rdBase, RDConfig
from rdkit.Chem import AllChem, rdShapeHelpers
from rdkit.Chem.FeatMaps import FeatMaps


### start function definitions #########################################

def log(*args, **kwargs):
    """Log output to STDERR
    """
    print(*args, file=sys.stderr, **kwargs)


# Setting up the features to use in FeatureMap
fdef = AllChem.BuildFeatureFactory(os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef'))

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
    #log("Calc:", str(fms[0].ScoreFeats(featLists[1])), "/ float(min(", str(fms[0].GetNumFeatures()), str(len(featLists[1])), "))")
    fm_score = fms[0].ScoreFeats(featLists[1]) / float(min(fms[0].GetNumFeatures(), len(featLists[1])))
    return fm_score

def get_SucosScore(ref_mol, query_mol, field_name):
    fm_score = get_FeatureMapScore(ref_mol, query_mol)
    #log("FeatureMapScore:", str(fm_score))
    protrude_dist = rdShapeHelpers.ShapeProtrudeDist(ref_mol, query_mol, allowReordering=False)
    #log("ProtrudeDistance:", str(protrude_dist))
    #log("Sucos calc: 0.5 *", str(fm_score), "+ 0.5 * (1.0 -", protrude_dist, ")")
    score = 0.5 * fm_score + 0.5 * (1.0 - protrude_dist)
    #log("SucosScore:", str(score))
    query_mol.SetDoubleProp(field_name, score)
    return score

def read_single_molecule(filename, index=1, format=None):
    """Read a single molecule as a RDKit Mol object. This can come from a file in molfile or SDF format.
    If SDF then you can also specify an index of the molecule that is read (default is the first)
    """
    mol = None
    if format == 'mol' or filename.lower().endswith('.mol') or filename.lower().endswith('.mol.gz'):
        file = open_file_for_reading(filename)
        mol = Chem.MolFromMolBlock(file.read())
        file.close()
    elif format == 'sdf' or filename.lower().endswith('.sdf') or filename.lower().endswith('.sdf.gz'):
        file = open_file_for_reading(filename)
        supplier = Chem.ForwardSDMolSupplier(file)
        for i in range(0,index):
            if supplier.atEnd():
                break
            mol = next(supplier)
        file.close()

    if not mol:
        raise ValueError("Unable to read molecule")

    return mol

def open_file_for_reading(filename):
    """Open the file gunzipping it if it ends with .gz."""
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'rb')
    else:
        return open(filename, 'rb')

def open_file_for_writing(filename):
    if filename.lower().endswith('.gz'):
        return gzip.open(filename, 'at')
    else:
        return open(filename, 'w+')

### start main execution #########################################

def main():

    parser = argparse.ArgumentParser(description='SuCOS with RDKit')
    parser.add_argument('-i', '--input', help='Input file in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('-o', '--output', help='Output file in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('-t', '--target', help='Molecule to compare against on Molfile (.mol) or SDF (.sdf) format')
    parser.add_argument('--target-format', help='Format for the target (mol or sdf)')
    parser.add_argument('--targetidx', help='Target molecule index in SD file if not the first', type=int, default=1)

    args = parser.parse_args()
    log("SuCOS Args: ", args)


    ref_mol = read_single_molecule(args.target, index=args.targetidx, format=args.target_format)
    log("Reference mol has", str(ref_mol.GetNumHeavyAtoms()), "heavy atoms")

    input_file = open_file_for_reading(args.input)
    suppl = Chem.ForwardSDMolSupplier(input_file)
    output_file = open_file_for_writing(args.output)
    writer = Chem.SDWriter(output_file)

    count = 0
    total = 0
    errors = 0
    for mol in suppl:
        count +=1
        if mol is None:
            continue
        #log("Mol has", str(mol.GetNumHeavyAtoms()), "heavy atoms")
        try:
            fm_score = get_SucosScore(ref_mol, mol, "SuCOS_Score")
            log("Score:", str(fm_score))
            writer.write(mol)
            total +=1
        except ValueError as e:
            errors +=1
            log("Molecule", count, "failed to score:", e.message)

    input_file.close()
    writer.flush()
    writer.close()
    output_file.close()

    log("Completed.", total, "processed, ", count, "succeeded, ", errors, "errors")

if __name__ == "__main__":
    main()
