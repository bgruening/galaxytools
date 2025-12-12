#!/usr/bin/env python

import sys
from pathlib import Path

import pandas as pd


def get_column_names(file_path, indices):
    """
    Get column names based on multiple indices from a tabular file.
    """
    try:
        file_ext = Path(file_path).suffix.lower()
        sep = "," if file_ext == ".csv" else "\t"

        if file_ext in [".csv", ".tsv", ".txt", ".tab", ".tabular"]:
            data = pd.read_csv(file_path, sep=sep)
        else:
            raise ValueError(f"Unsupported file extension: {file_ext}")

    except Exception as e:
        raise ValueError(f"Error loading data from {file_path}: {e}") from e

    column_names = []

    for index in indices:
        zero_based_index = index - 1  # Convert to zero-based index
        if zero_based_index < 0 or zero_based_index >= len(data.columns):
            print(
                f"Error: Index {index} is out of range. File has {len(data.columns)} columns (1-{len(data.columns)})."
            )
            return None
        column_names.append(data.columns[zero_based_index].strip())

    return column_names


if __name__ == "__main__":

    file_path = sys.argv[1]
    indices_str = sys.argv[2]

    # Parse comma-separated indices
    try:
        indices = [int(i.strip()) for i in indices_str.split(",")]
    except ValueError:
        print("Error: All indices must be integers.")
        sys.exit(1)

    result = get_column_names(file_path, indices)
    if result is not None:
        print(",".join(result))
