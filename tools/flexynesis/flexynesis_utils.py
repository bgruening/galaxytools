#!/usr/bin/env python

import os
import argparse
import pandas as pd
import sys
from pathlib import Path


def make_data_dict(clin_path, omics_paths):
    """Read clinical and omics data files into a dictionary."""
    data = {}

    # Read clinical data
    print(f"Reading clinical data from {clin_path}")
    try:
        file_ext = Path(clin_path).suffix.lower()
        if file_ext == '.csv':
            clin = pd.read_csv(clin_path, index_col=0)
        elif file_ext in ['.tsv', '.txt', '.tab', '.tabular']:
            clin = pd.read_csv(clin_path, sep='\t', index_col=0)
        else:
            raise ValueError(f"Unsupported clinical file format: {file_ext}")

        if clin.empty:
            raise ValueError(f"Clinical file {clin_path} is empty")
        data['clin'] = clin
        print(f"Loaded clinical data: {clin.shape[0]} samples, {clin.shape[1]} features")
    except Exception as e:
        raise ValueError(f"Error reading clinical file {clin_path}: {e}")

    # Read omics data
    print(f"Reading omics data from {', '.join(omics_paths)}")
    for path in omics_paths:
        try:
            name = os.path.splitext(os.path.basename(path))[0]
            file_ext = Path(path).suffix.lower()
            if file_ext == '.csv':
                df = pd.read_csv(path, index_col=0)
            elif file_ext in ['.tsv', '.txt', '.tab', '.tabular']:
                df = pd.read_csv(path, sep='\t', index_col=0)
            else:
                raise ValueError(f"Unsupported omics file format: {file_ext}")
            if df.empty:
                print(f"Warning: Omics file {path} is empty, skipping")
                continue
            data[name] = df
            print(f"Loaded {name}: {df.shape[0]} features, {df.shape[1]} samples")
        except Exception as e:
            print(f"Warning: Error reading omics file {path}: {e}")
            continue

    if len(data) == 1:  # Only clinical data loaded
        raise ValueError("No omics data was successfully loaded")

    return data


def validate_data_consistency(data):
    """Validate that clinical and omics data have consistent samples."""
    clin_samples = set(data['clin'].index)

    for name, df in data.items():
        if name == 'clin':
            continue

        omics_samples = set(df.columns)

        # Check for sample overlap
        common_samples = clin_samples.intersection(omics_samples)
        if len(common_samples) == 0:
            raise ValueError(f"No common samples between clinical data and {name}")

        missing_in_omics = clin_samples - omics_samples
        missing_in_clin = omics_samples - clin_samples

        if missing_in_omics:
            print(f"Warning: {len(missing_in_omics)} clinical samples not found in {name}")
        if missing_in_clin:
            print(f"Warning: {len(missing_in_clin)} samples in {name} not found in clinical data")

    return True


def split_and_save_data(data, ratio=0.7, output_dir='.'):
    """Split data into train/test sets and save to files."""
    # Validate data consistency first
    validate_data_consistency(data)

    samples = data['clin'].index.tolist()

    train_samples = list(pd.Series(samples).sample(frac=ratio, random_state=42))
    test_samples = list(set(samples) - set(train_samples))

    train_data = {}
    test_data = {}

    for key, df in data.items():
        try:
            if key == 'clin':
                train_data[key] = df.loc[df.index.intersection(train_samples)]
                test_data[key] = df.loc[df.index.intersection(test_samples)]
            else:
                train_data[key] = df.loc[:, df.columns.intersection(train_samples)]
                test_data[key] = df.loc[:, df.columns.intersection(test_samples)]
        except Exception as e:
            print(f"Error splitting data {key}: {e}")
            continue

    # Create output directories
    os.makedirs(os.path.join(output_dir, 'train'), exist_ok=True)
    os.makedirs(os.path.join(output_dir, 'test'), exist_ok=True)

    # Save train and test data
    for key in data.keys():
        try:
            train_data[key].to_csv(os.path.join(output_dir, 'train', f'{key}.csv'))
            test_data[key].to_csv(os.path.join(output_dir, 'test', f'{key}.csv'))
        except Exception as e:
            print(f"Error saving {key}: {e}")
            continue


def main():
    parser = argparse.ArgumentParser(description='Flexynesis splitting pipeline')
    parser.add_argument('--clin', required=True,
                        help='Path to clinical data CSV file (samples in rows)')
    parser.add_argument('--omics', required=True,
                        help='Comma-separated list of omics CSV files (samples in columns)')
    parser.add_argument('--split', type=float, default=0.7,
                        help='Train split ratio (default: 0.7)')
    parser.add_argument('--out', default='.',
                        help='Output directory (default: current directory)')

    args = parser.parse_args()

    try:
        # Validate inputs
        if not os.path.isfile(args.clin):
            raise FileNotFoundError(f"Clinical file not found: {args.clin}")

        # Parse omics files
        omics_files = [f.strip() for f in args.omics.split(',') if f.strip()]
        if not omics_files:
            raise ValueError("At least one omics file must be provided")

        # Check omics files exist
        for f in omics_files:
            if not os.path.isfile(f):
                raise FileNotFoundError(f"Omics file not found: {f}")

        # Validate split ratio
        if not 0 < args.split < 1:
            raise ValueError(f"Split ratio must be between 0 and 1, got {args.split}")

        # Create output directory if it doesn't exist
        if not os.path.exists(args.out):
            os.makedirs(args.out)

        data = make_data_dict(args.clin, omics_files)
        split_and_save_data(data, ratio=args.split, output_dir=args.out)

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
