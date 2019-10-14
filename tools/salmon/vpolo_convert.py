#!/usr/bin/env python

from vpolo.alevin import parser as par
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("--mtx", "-m", action="store_true", help="--dumpMtx flag set")
parser.add_argument("--umi", "-u", action="store_true", help="--dumpUmiGraph flag set")
args = parser.parse_args()

if args.mtx:
    alevin_df = par.read_quants_bin("output")
    with open("quants_mat.tsv", "w") as f:
        f.write(alevin_df.to_csv(sep="\t"))

if args.umi:
    os.mkdir("umiout")
    par.read_umi_graph("output", "umiout")
