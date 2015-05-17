#!/usr/bin/env python

"""
Takes a ciros data file ...

chromosom <tab> start <tab> end <tab> additinal information

and reformate it. So that the start column is really smaller than the end column.

Example call:
./reformatCircosDataFile_start_smaller_than_end.py input_file.dat output_file.dat 2 3

"""

import sys

in_file = sys.argv[1]
out_file = open(sys.argv[2], 'w')

start_col = int(sys.argv[3])
end_col = int(sys.argv[4])


for line in open( in_file ):
    cols = line.strip().split()
    if not line.strip():
        continue
    start = cols[start_col - 1]
    end = cols[ end_col - 1 ]
    if int(start) > int(end):
        cols[start_col - 1] = end
        cols[ end_col - 1 ] = start
        out_file.write( '\t'.join(cols) + '\n' )
    else:
        out_file.write( line )

out_file.close()
