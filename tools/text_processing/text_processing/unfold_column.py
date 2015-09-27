#!/usr/bin/env python

import sys

out = open(sys.argv[4], 'w+')

sep = sys.argv[3]
# un-sanitize Galaxy inputs
if sep == 'X':
    sep = ';'

with open(sys.argv[1]) as handle:
    for line in handle:
        cols = line.split('\t')
        unfolding_column = int(sys.argv[2]) - 1
        column_content = cols[ unfolding_column ]
        for elem in column_content.split( sep ):
            out.write( '\t'.join( cols[:unfolding_column] + [elem] + cols[unfolding_column+1:]) )
out.close()
