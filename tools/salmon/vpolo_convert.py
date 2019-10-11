#!/usr/bin/env python

from vpolo.alevin import parser
#import argparse

#parser = argparse.ArgumentParser()
#parser.add_argument("--mtx", "-m", action="store_true", help="--dumpMtx flag set")
#parser.add_argument("--umi", "-u", action="store_true", help="--dumpUmiGraph flag set")
#args = parser.parse_args()

if args.mtx:
    alevin_df = parser.read_quants_bin("output")
    with open("quants_mat.tsv", "w") as f:
        f.write(alevin_df.to_csv(sep="\t"))

#if args.umi:
#    parser.read_umi_graph("output", "umigraph")
