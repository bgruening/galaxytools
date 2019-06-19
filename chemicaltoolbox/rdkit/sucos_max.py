#!/usr/bin/env python
"""
Determination of the best SuCOS score for one set of molecules against a second set.
This is a quite specialised function that is designed to take a set of potential follow up
ligands and compare them to a set of clustered fragment hits to help identify which follow up
ligands best map to the binding space of the hits.

The clustering if the fragment hits is expected to be performed with the sucos_cluster.py module
and will generate a set of SD files, one for each cluster of hits (presumably corresponding to a
binding pocket in the protein target).

Each molecule in the input ligands is then compared (using SuCOS) to each hit in the clusters to
identify the hit with the best SuCOS score. The output is a SD file with each of the ligands, with
these additional fields for each molecule:
Max_SuCOS_Score - the best score
Max_SuCOS_FeatureMap_Score - the feature map score for the hit that has the best SuCOS score
Max_SuCOS_Tanimoto_Score - the Tanimoto score for the hit that has the best SuCOS score
Max_SuCOS_Cluster - the name of the cluster SD file that contains the best hit
Max_SuCOS_Index - the index of the best hit in the SD file

If a molecule has no alignment to any of the clustered hits (a max score of zero) then it is not
included in the results.


SuCOS is the work of Susan Leung.
Bitbucket: https://bitbucket.org/Susanhleung/sucos/
Publication: https://doi.org/10.26434/chemrxiv.8100203.v1
"""

import sucos, utils
import argparse, gzip
from rdkit import Chem


def process(inputfilename, clusterfilenames, outputfilename):

    all_clusters = {}
    for filename in clusterfilenames:
        cluster = []
        cluster_file = utils.open_file_for_reading(filename)
        suppl = Chem.ForwardSDMolSupplier(cluster_file)
        for mol in suppl:
            if not mol:
                continue
            features = sucos.getRawFeatures(mol)
            cluster.append((mol, features))
        cluster_file.close()
        all_clusters[filename] = cluster

    input_file = utils.open_file_for_reading(inputfilename)
    suppl = Chem.ForwardSDMolSupplier(input_file)
    output_file = utils.open_file_for_writing(outputfilename)
    writer = Chem.SDWriter(output_file)

    comparisons = 0
    mol_num = 0

    for mol in suppl:
        mol_num += 1
        max_sucos_score = 0
        cluster_name = None
        cluster_index = 0
        for clusterfilename in all_clusters:
            cluster = all_clusters[clusterfilename]
            index = 0
            for entry in cluster:
                hit = entry[0]
                features = entry[1]
                index += 1
                comparisons += 1
                sucos_score, fm_score, tani_score = sucos.get_SucosScore(hit, mol, ref_features=features)
                if sucos_score > max_sucos_score:
                    max_sucos_score = sucos_score
                    max_fm_score = fm_score
                    max_tanimoto_score = tani_score
                    cluster_name = clusterfilename
                    cluster_index = index

        utils.log("Max SuCOS:", max_sucos_score, "File:", cluster_name, "Index:", cluster_index)
        if max_sucos_score > 0:
            mol.SetDoubleProp("Max_SuCOS_Score", max_sucos_score)
            mol.SetDoubleProp("Max_SuCOS_FeatureMap_Score", max_fm_score)
            mol.SetDoubleProp("Max_SuCOS_Tanimoto_Score", max_tanimoto_score)
            mol.SetProp("Max_SuCOS_Cluster", cluster_name)
            mol.SetIntProp("Max_SuCOS_Index", cluster_index)
            writer.write(mol)
        else:
            utils.log("Molecule", mol_num, "did not overlay. Omitting from results")


    input_file.close()
    writer.flush()
    writer.close()
    output_file.close()

    utils.log("Completed", comparisons, "comparisons")


### start main execution #########################################

def main():
    parser = argparse.ArgumentParser(description='Max SuCOS scores with RDKit')
    parser.add_argument('-i', '--input', help='Input file to score in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('-o', '--output', help='Output file in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('clusters', nargs='*', help="One of more SDF files with the clustered hits")

    args = parser.parse_args()
    utils.log("Max SuCOS Args: ", args)

    process(args.input, args.clusters, args.output)


if __name__ == "__main__":
    main()