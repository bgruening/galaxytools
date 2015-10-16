#!/usr/bin/env python 
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Tabular Input File Name')
parser.add_argument('-o','--output', help='Tabular Output File')
parser.add_argument('-c', '--columns', nargs='+', help='Column Headers to Sort By')
args=parser.parse_args()

cols = args.columns
table = pd.read_csv(args.input, sep='\t')
blist = list(table.columns)
for token in cols:
    blist.remove(token)
sorted_table = table[args.columns + blist]
# write without index, seperated by tabs
sorted_table.to_csv(args.output, sep='\t', index=False)
