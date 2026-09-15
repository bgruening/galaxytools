#!/usr/bin/env python

import sys
import csv
import numpy as np

import Bio.Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor

def usage():
  print("usage: mentalist_tree <input.tsv>\n")

def process_input_matrix(input_matrix):
    """ Converts an array-of-arrays containting sample IDs and distances
        into a BioPython DistanceMatrix object
    """
    input_matrix.pop(0)
    sample_names = [row[0] for row in input_matrix]
    for row in input_matrix:
        row.pop(0)
    distance_matrix = []
    for input_matrix_row in input_matrix:
        distance_matrix.append([int(i) for i in input_matrix_row])
    """ np.tril() converts a matrix like this: [[0 1 2]
                                                [1 0 1]
                                                [2 1 0]]
        ...into this: [[0 0 0]
                       [1 0 0]
                       [2 1 0]]
        ...but what we need to pass to DistanceMatrix() is this: [[0]
                                                                  [1 0]
                                                                  [2 1 0]]
        ...so that's what the (somewhat cryptic) code below does.
    """
    distance_matrix = np.tril(np.array(distance_matrix))
    num_rows = distance_matrix.shape[0]
    """ masking the distance matrix with tril_indices gives a linearized
        distance matrix [0 1 0 2 1 0] that we need to re-construct into [[0], [1, 0], [2, 1, 0]]
    """
    lower_triangular_idx_mask = np.tril_indices(num_rows)
    linear_distance_matrix = distance_matrix[lower_triangular_idx_mask]
    distance_matrix = []
    min = 0
    max = 1
    for i in range(num_rows):
        distance_matrix.append(linear_distance_matrix[min:max].tolist())
        min = max
        max = max + (i + 2) 
    distance_matrix = DistanceMatrix(names=sample_names, matrix=distance_matrix)
    return distance_matrix

def main():
    if len(sys.argv) < 2:
        usage()
        sys.exit(1)

    input_file = sys.argv[1]
    reader = csv.reader(open(input_file, "r"), delimiter="\t")
    input_matrix = list(reader)
    # Don't build a tree with fewer than 3 samples, just produce an empty file
    if len(input_matrix) < 4:
      print('();')
      sys.exit(0)
    distance_matrix = process_input_matrix(input_matrix)
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    Bio.Phylo.write(tree, sys.stdout, 'newick')

if __name__ == '__main__':
    main()
