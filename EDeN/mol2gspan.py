#!/usr/bin/env python

import os, sys
import argparse

def main(args ):

    begin = False
    iid = 0
    graph_counter = 1

    for line in args.infile:
        if line.rstrip():
            if line.strip().endswith('END'):
                begin = False
            elif line.strip() == '$$$$':
                graph_counter += 1
                iid = 0
            else:
                # found header line, like:  21 21  0  0  0  0  0  0  0  0999 V2000
                if len(line.split()) >= 5 and line.split()[-1] == 'V2000':
                    args.outfile.write('t # id %s\n' % graph_counter)
                    begin=True
                    continue
                # connection or coordinate/atom table
                if len(line.split()) >= 4 and begin:
                    # coordinate/atom table
                    if not line.startswith('M'):
                        if line.split()[3].isalpha() or line.split()[3] == '*':
                            args.outfile.write( 'v %s %s \n' % (iid, line.split()[3]) )
                            iid += 1
                        else:
                            #connection table
                            id, node, edge, trash = line.split(None, 3)
                            args.outfile.write( 'e %s %s %s\n' % ( int(id) - 1 , int(node) -1, edge ) )

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', nargs='?', type=argparse.FileType('r'),
        default=sys.stdin, help="Specify one or more input files")
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        default=sys.stdout, help="Specify one output file")
    args = parser.parse_args()
    main( args )
