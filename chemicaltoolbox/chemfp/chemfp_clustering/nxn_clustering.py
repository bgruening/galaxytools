#!/usr/bin/env python
"""
    Modified version of code examples from the chemfp project.
    http://code.google.com/p/chem-fingerprints/
    Thanks to Andrew Dalke of Andrew Dalke Scientific!
"""
import matplotlib
matplotlib.use('Agg')
import argparse
import os
import chemfp
import scipy.cluster.hierarchy as hcluster
import pylab
import numpy


def distance_matrix(arena, tanimoto_threshold = 0.0):
    n = len(arena)
    # Start off a similarity matrix with 1.0s along the diagonal
    try:
        similarities = numpy.identity(n, "d")
    except:
        raise Exception('Input dataset is to large!')
    chemfp.set_num_threads( args.processors )

    ## Compute the full similarity matrix.
    # The implementation computes the upper-triangle then copies
    # the upper-triangle into lower-triangle. It does not include
    # terms for the diagonal.
    results = chemfp.search.threshold_tanimoto_search_symmetric(arena, threshold=tanimoto_threshold)

    # Copy the results into the NumPy array.
    for row_index, row in enumerate(results.iter_indices_and_scores()):
        for target_index, target_score in row:
            similarities[row_index, target_index] = target_score

    # Return the distance matrix using the similarity matrix
    return 1.0 - similarities



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="""NxN clustering for fps files.
For more details please see the chemfp documentation:
https://chemfp.readthedocs.org
""")

    parser.add_argument("-i", "--input", dest="input_path",
                    required=True,
                    help="Path to the input file.")

    parser.add_argument("-o", "--output", dest="output_path",
                    help="Path to the output file.")

    parser.add_argument("-t", "--threshold", dest="tanimoto_threshold", 
                    type=float, default=0.0,
                    help="Tanimoto threshold [0.0]")

    parser.add_argument("--oformat", default='png', help="Output format (png, svg).")

    parser.add_argument('-p', '--processors', type=int, 
        default=4)

    args = parser.parse_args()


    arena = chemfp.load_fingerprints( args.input_path )
    distances  = distance_matrix( arena, args.tanimoto_threshold )
    linkage = hcluster.linkage( distances, method="single", metric="euclidean" )

    hcluster.dendrogram(linkage, labels=arena.ids)

    pylab.savefig( args.output_path, format=args.oformat )



