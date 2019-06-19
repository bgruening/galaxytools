#!/usr/bin/env python
"""
Cluster a set of molecules based on their 3D overlays as determined by the SuCOS score.

This will generate a set of SD files, one for each cluster of molecules (presumably corresponding to a
binding pocket in the protein target).


SuCOS is the work of Susan Leung.
Bitbucket: https://bitbucket.org/Susanhleung/sucos/
Publication: https://doi.org/10.26434/chemrxiv.8100203.v1
"""

import sucos, utils
import argparse, gzip
from rdkit import Chem
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import linkage, fcluster

### start main execution #########################################


def calc_distance_matrix(mols):
    """
    Calculate a full distance matrix for the given molecules. Identical molecules get a score of 0.0 with the maximum
    distance possible being 1.0.
    :param mols: A list of molecules. It must be possible to iterate through this list multiple times
    :return: A NxN 2D array of distance scores, with N being the number of molecules in the input
    """

    # TODO - do we need to calculate both sides of the matrix? Tanimoto is supposed to be a symmetric distance measure,
    #  but the matrix that is generated does not seem to be symmetric.

    mol_fm_tuples = []
    for mol in mols:
        features = sucos.getRawFeatures(mol)
        mol_fm_tuples.append((mol, features))

    matrix = []
    for tuple1 in mol_fm_tuples:
        tmp = []
        for tuple2 in mol_fm_tuples:
            if tuple1[0] == tuple2[0]:
                tmp.append(0.0)
            else:
                #utils.log("Calculating SuCOS between", mol1, mol2)
                sucos_score, fm_score, tani_score = sucos.get_SucosScore(tuple1[0], tuple2[0],
                    tani=True, ref_features=tuple1[1], query_features=tuple2[1])
                tmp.append(1.0 - sucos_score)
        matrix.append(tmp)


    return matrix


def cluster(matrix, threshold=0.8):
    """
    Cluster the supplied distance matrix returning an array of clusters.
    :param matrix: the distance matrix, as calculated with the calc_distance_matrix function.
    :param threshold: The clustering cuttoff. The default of 0.8 is a reasonable value to use.
    :return: An array of clusters, each cluster being an array of the indices from the matrix.
    """

    indexes = [x for x in range(0, len(matrix))]
    cols = [x for x in range(0, len(matrix[0]))]
    #utils.log("indexes", indexes)
    #utils.log("cols", cols)
    df = pd.DataFrame(matrix, columns=cols, index=indexes)
    utils.log("DataFrame:", df.shape)
    #utils.log(df)
    indices = np.triu_indices(df.shape[0], k=1)
    #utils.log("Indices:", indices)
    t = np.array(df)[indices]
    Z = linkage(t, 'average')
    lig_clusters = []
    cluster_arr = fcluster(Z, t=threshold, criterion='distance')
    for i in range(np.amax(cluster_arr)):
        clus = df.columns[np.argwhere(cluster_arr==i+1)]
        lig_clusters.append([x[0] for x in clus.tolist()])

    utils.log("Clusters", lig_clusters)
    return lig_clusters

def write_clusters_to_sdfs(mols, clusters, basename, gzip=False):
    """
    Write the molecules to SDF files, 1 file for each cluster.
    :param mols The molecules to write:
    :param clusters The clusters, as returned by the cluster function:
    :param basename The basename for the file name. e.g. if basename is 'output' then files like
    output1.sdf, output2.sdf will be written:
    :param gzip Whether to gzip the output
    :return:
    """

    i = 0
    for cluster in clusters:
        i += 1
        filename = basename + str(i) + ".sdf"
        if gzip:
            filename += ".gz"
        utils.log("Writing ", len(cluster), "molecules in cluster", i, "to file", filename)
        output_file = utils.open_file_for_writing(filename)
        writer = Chem.SDWriter(output_file)
        for index in cluster:
            mol = mols[index]
            writer.write(mol)
        writer.flush()
        writer.close()
        output_file.close()



def main():
    parser = argparse.ArgumentParser(description='Clustering with SuCOS and RDKit')
    parser.add_argument('-i', '--input', help='Input file in SDF format. Can be gzipped (*.gz).')
    parser.add_argument('-o', '--output', default="cluster", help="Base name for output files in SDF format. " +
                                               "e.g. if value is 'output' then files like output1.sdf, output2.sdf will be created")
    parser.add_argument('--gzip', action='store_true', help='Gzip the outputs generating files like output1.sdf.gz, output2.sdf.gz')
    parser.add_argument('-t', '--threshold', type=float, default=0.8, help='Clustering threshold')

    args = parser.parse_args()
    utils.log("SuCOS Cluster Args: ", args)

    input_file = utils.open_file_for_reading(args.input)
    suppl = Chem.ForwardSDMolSupplier(input_file)
    mols = list(suppl)
    matrix = calc_distance_matrix(mols)
    clusters = cluster(matrix, threshold=args.threshold)
    write_clusters_to_sdfs(mols, clusters, args.output, gzip=args.gzip)


if __name__ == "__main__":
    main()