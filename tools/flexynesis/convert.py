#!/usr/bin/env python

import sys

import pandas as pd


def tabular_to_csv(tabular_file, csv_file):
    """Convert tabular (TSV) to CSV"""
    data = pd.read_csv(tabular_file, sep="\t")
    if data.columns[0] == '' or data.columns[0].startswith('Unnamed:'):
        data.columns = ['ID'] + list(data.columns[1:])
    data.to_csv(csv_file, index=False)


def csv_to_tabular(csv_file, tabular_file):
    """Convert CSV to tabular (TSV)"""
    data = pd.read_csv(csv_file)
    if data.columns[0] == '' or data.columns[0].startswith('Unnamed:'):
        data.columns = ['ID'] + list(data.columns[1:])
    data.to_csv(tabular_file, sep="\t", index=False)


if __name__ == "__main__":
    input_file = sys.argv[1]
    output_file = sys.argv[2]

    if input_file.endswith('.csv'):
        csv_to_tabular(input_file, output_file)
    else:
        tabular_to_csv(input_file, output_file)
