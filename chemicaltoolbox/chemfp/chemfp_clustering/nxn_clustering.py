#!/usr/bin/env python
"""
    Modified version of code examples from the chemfp project.
    http://code.google.com/p/chem-fingerprints/
    Thanks to Andrew Dalke of Andrew Dalke Scientific!
"""
import matplotlib
matplotlib.use('Agg')
import sys
import os
import chemfp
import scipy.cluster.hierarchy as hcluster
import pylab
import numpy


def distance_matrix(arena,t):
    n = len(arena)
    # The Tanimoto search computes all of the scores when threshold=0.0.
    # The SearchResult contains sparse data, so I set all values
    # now to 1.0 so you can experiment with higher thresholds.
    distances = numpy.ones((n, n), numpy.float64)

    # Keep track of where the query subarena is in the query
    query_row = 0

    for query_arena in arena.iter_arenas():
        results = arena.threshold_tanimoto_search_arena(query_arena, threshold=t)  
    for q_i, hits in enumerate(results.iter_indices_and_scores()):
            query_idx = query_row + q_i
            for target_idx, score in hits:
                distances[query_idx, target_idx] = 1.0 - score
        query_row += len(query_arena)

    return distances

dataset = chemfp.load_fingerprints( sys.argv[1] )
distances  = distance_matrix( dataset,float( sys.argv[2] ) )
linkage = hcluster.linkage( distances, method="single", metric="euclidean" )

# Plot using matplotlib, which you must have installed
hcluster.dendrogram(linkage, labels=dataset.ids)

pylab.savefig( sys.argv[3], format='svg' )







