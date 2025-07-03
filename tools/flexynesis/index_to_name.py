#!/usr/bin/env python

import sys

import pandas as pd

def get_column_name(file_path, index):
    """
    Get the column name based on an index from a tabular file.
    """
    data = pd.read_csv(file_path, sep=",")
    index = index-1  # Convert to zero-based index
    if index < 0 or index >= len(data.columns):
        print(f"Error: Index {index+1} is out of range. File has {len(data.columns)} columns (1-{len(data.columns)}).")
        return None
    return data.columns[index].strip()

if __name__ == "__main__":

    file_path = sys.argv[1]
    index = int(sys.argv[2])

    print(get_column_name(file_path, index))
