#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Tabular Input File Name')
parser.add_argument('-o','--output', help='Tabular Output File')
parser.add_argument(
    '-c', '--columns', nargs='+', help='Column Headers to Sort By'
)
parser.add_argument(
    '-d', '--discard', action='store_true',
    help='Discard remaining columns'
)

args=parser.parse_args()

with open(args.input) as data:
    hdr = next(data)
    columns = hdr.rstrip('\n').split('\t')
    idx = [columns.index(name) for name in args.columns]
    if not args.discard:
        idx += [i for i in range(len(columns)) if i not in idx]
    rearranged_cols = [columns[i] for i in idx]
    with open(args.output, 'w') as out:
        out.write('\t'.join(rearranged_cols) + '\n')
        for line in data:
            columns = line.rstrip('\n').split('\t')
            rearranged_cols = [columns[i] for i in idx]
            out.write('\t'.join(rearranged_cols) + '\n')
