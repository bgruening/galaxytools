#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import re
from pathlib import Path

taxo_level = {
    'd': 'domain',
    'k': 'kingdom',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
    't': 'strains'}


def format_for_krona(metaphlan_output_fp, krona_out_fp):
    '''
    Split default MetaPhlAn into a report for each taxonomic levKRONAel

    :param metaphlan_output_fp: Path default MetaPhlAn output
    :param krona_out: Path to output file for Krona
    '''
    re_replace = re.compile(r"\w__")
    re_bar = re.compile(r"\|")
    re_underscore = re.compile(r"_")

    with open(metaphlan_output_fp, 'r') as metaphlan_output_f:
        with open(krona_out_fp, 'w') as krona_out_f:
            for line in metaphlan_output_f.readlines():
                if "s__" in line:
                    x = line.rstrip().split('\t')
                    x.pop()
                    lineage = re.sub(re_bar, '', x[0])
                    lineage = re.sub(re_replace, '\t', lineage)
                    lineage = re.sub(re_underscore, ' ', lineage)
                    krona_out_f.write("%s\t%s\n" % (x[-1], lineage))


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Format MetaPhlAn output')
    subparsers = parser.add_subparsers(dest='function')
    # split_levels
    split_levels_parser = subparsers.add_parser('split_levels', help='Split default MetaPhlAn into a report for each taxonomic level')
    split_levels_parser.add_argument('--metaphlan_output', help="Path to default MetaPhlAn output")
    split_levels_parser.add_argument('--outdir', help="Path to output directory")
    split_levels_parser.add_argument('--legacy-output', dest='legacy_output', action='store_true', help="Old MetaPhlAn2 two columns output")
    split_levels_parser.set_defaults(legacy_output=False)
    # format_for_krona
    format_for_krona_parser = subparsers.add_parser('format_for_krona', help='Split default MetaPhlAn into a report for each taxonomic level')
    format_for_krona_parser.add_argument('--metaphlan_output', nargs='*', help="Path to default MetaPhlAn output")
    format_for_krona_parser.add_argument('--krona_output', help="Path to Krona output directory")

    args = parser.parse_args()
    krona_out = args.krona_output
    for i in args.metaphlan_output:
        krona_outfile = i + "_" + krona_out
        format_for_krona(
            Path(i),
            Path(krona_outfile))
