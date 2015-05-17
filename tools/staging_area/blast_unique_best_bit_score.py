#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys
import argparse


def stop_err( msg ):
    sys.stderr.write("%s\n" % msg)
    sys.exit(1)

def get_col_index(col_str):
    if isinstance(col_str, str) and col_str[0]=="c":
        col_str = col_str[1:]
    return int(col_str)-1

def main( blast_tabular_file, outfile, score_column = '12', sort_order = 'high'):

    want_highest = want_lowest = False
    if sort_order == "high":
        want_highest = True
    elif sort_order == "low":
        want_lowest = True
    else:
        stop_err("Sort order argument should be high or low")

    score_column = get_col_index(score_column)

    ofile = open(outfile, 'w+')
    best_hits = {}
    for line in open(blast_tabular_file):
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        cols = line.split('\t')
        query_id = cols[0]

        score = float(cols[ score_column ])
        if not best_hits.has_key(query_id):
            best_hits[query_id] = (score, line)
        elif want_highest and best_hits[query_id][0] < score:
            best_hits[query_id] = (score, line)
        elif want_lowest and best_hits[query_id][0] > score:
            best_hits[query_id] = (score, line)

    for key in sorted(best_hits.keys()):
        ofile.write( best_hits[key][1] + '\n' )
    ofile.close()





if __name__ == "__main__" :
    main(sys.argv[1], sys.argv[2], int(sys.argv[3]), sys.argv[4])
