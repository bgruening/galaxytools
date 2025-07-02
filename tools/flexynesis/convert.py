#!/usr/bin/env python

import sys

import pandas as pd


def tabular_to_csv(tabular_file, csv_file):
    """Convert tabular (TSV) to CSV"""
    data = pd.read_csv(tabular_file, sep="\t")
    data.to_csv(csv_file, index=False)


if __name__ == "__main__":
    tabular_to_csv(sys.argv[1], sys.argv[2])
