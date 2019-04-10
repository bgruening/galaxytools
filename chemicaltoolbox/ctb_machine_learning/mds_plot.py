#!/usr/bin/env python

import argparse
import os
import sklearn.manifold
import numpy
import math
import pylab

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""2D multidimensional scaling of NxN matrices with scatter plot"""
        )

    parser.add_argument("-i", "--input", dest="sm",
                    required=True,
                    help="Path to the input file.")
    parser.add_argument("--oformat", default='png', help="Output format (png, svg)")
    parser.add_argument("-o", "--output", dest="output_path",
                    help="Path to the output file.")

    args = parser.parse_args()
    mds = sklearn.manifold.MDS( n_components=2, max_iter=300, eps=1e-6, dissimilarity='precomputed' )
    data = numpy.loadtxt(args.sm)
    d = len(data)
    sm = numpy.reshape( data, ( d,d ))
    pos = mds.fit( sm ).embedding_
    pylab.scatter( pos[:,0],pos[:,1] )
    pylab.axis('off')
    pylab.savefig( args.output_path, format=args.oformat )
