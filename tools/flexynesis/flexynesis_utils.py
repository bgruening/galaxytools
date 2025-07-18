#!/usr/bin/env python

import argparse
import os
import sys
from pathlib import Path

import pandas as pd


def read_data(data_input, index=False):
    """Load CSV or TSV data file."""
    try:
        file_ext = Path(data_input).suffix.lower()
        sep = ',' if file_ext == '.csv' else '\t'
        index_col = 0 if index else None

        if file_ext in ['.csv', '.tsv', '.txt', '.tab', '.tabular']:
            return pd.read_csv(data_input, sep=sep, index_col=index_col)
        else:
            raise ValueError(f"Unsupported file extension: {file_ext}")
    except Exception as e:
        raise ValueError(f"Error loading data from {data_input}: {e}") from e


def binarize_mutations(df, gene_idx=1, sample_idx=2):
    """
    Binarize mutation data by creating a matrix of gene x sample with 1/0 values.
    """
    # galaxy index is 1-based, convert to zero-based
    gene_idx -= 1
    sample_idx -= 1
    # check idx
    if gene_idx >= len(df.columns) or sample_idx >= len(df.columns):
        raise ValueError(f"Column indices out of bounds. DataFrame has {len(df.columns)} columns, "
                         f"but requested indices are {gene_idx} and {sample_idx}")
    if gene_idx == sample_idx:
        raise ValueError("Gene and sample column indices must be different")

    # Get column names by index
    gene_col = df.columns[gene_idx]
    print(f"Using gene column: {gene_col} (index {gene_idx})")
    sample_col = df.columns[sample_idx]
    print(f"Using sample column: {sample_col} (index {sample_idx})")

    # Check if columns contain data
    if df[gene_col].isna().all():
        raise ValueError(f"Gene column (index {gene_idx}) contains only NaN values.")
    if df[sample_col].isna().all():
        raise ValueError(f"Sample column (index {sample_idx}) contains only NaN values.")

    # Group by gene and sample, count mutations
    mutation_counts = df.groupby([gene_col, sample_col]).size().reset_index(name='count')

    # Create pivot table
    mutation_matrix = mutation_counts.pivot(index=gene_col, columns=sample_col, values='count').fillna(0)

    # Binarize: convert any count > 0 to 1
    mutation_matrix[mutation_matrix > 0] = 1

    return mutation_matrix


def make_data_dict(clin_path, omics_paths):
    """Read clinical and omics data files into a dictionary."""
    data = {}

    # Read clinical data
    print(f"Reading clinical data from {clin_path}")
    try:
        clin = read_data(clin_path, index=True)

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
            df = read_data(path, index=True)
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
            train_data[key].to_csv(os.path.join(output_dir, 'train', f'{key}.tabular'), sep='\t')
            test_data[key].to_csv(os.path.join(output_dir, 'test', f'{key}.tabular'), sep='\t')
        except Exception as e:
            print(f"Error saving {key}: {e}")
            continue


def validate_survival(df, column_name):
    """Validate that survival column in the DataFrame contains numeric integers."""
    if column_name not in df.columns:
        raise ValueError(f"Column '{column_name}' not found in DataFrame.")

    try:
        numeric_col = pd.to_numeric(df[column_name], errors='raise')
    except Exception as e:
        raise ValueError(f"Non-numeric values found in column '{column_name}': {e}")

    if not (numeric_col.dropna() == numeric_col.dropna().astype(int)).all():
        raise ValueError(f"Column '{column_name}' contains non-integer numeric values.")

    print("All values are numeric integers.")


def main():
    parser = argparse.ArgumentParser(description='Flexynesis extra utilities')

    parser.add_argument("--util", type=str, required=True,
                        choices=['split', 'binarize', 'validate_survival'],
                        help="Utility function: 'split' for spiting data to train and test, 'binarize' for creating a binarized matrix from a mutation data, 'validate_survival' for validating survival data.")

    # Arguments for split (clin also for validate_survival)
    parser.add_argument('--clin', required=False,
                        help='Path to clinical data CSV file (samples in rows)')
    parser.add_argument('--omics', required=False,
                        help='Comma-separated list of omics CSV files (samples in columns)')
    parser.add_argument('--split', type=float, default=0.7,
                        help='Train split ratio (default: 0.7)')

    # Arguments for binarize
    parser.add_argument('--mutations', type=str, required=False,
                        help='Path to mutation data CSV file (samples in rows, genes in columns)')
    parser.add_argument('--gene_idx', type=int, default=0,
                        help='Column index for genes in mutation data (default: 0)')
    parser.add_argument('--sample_idx', type=int, default=1,
                        help='Column index for samples in mutation data (default: 1)')

    # Arguments for validate_survival
    parser.add_argument('--surv_event_var', type=str, required=False,
                        help='Column name for survival event variable (e.g., death)')

    # common arguments (binarize and split)
    parser.add_argument('--out', default='.',
                        help='Output directory (default: current directory)')

    args = parser.parse_args()

    try:
        # validate utility function
        if not args.util:
            raise ValueError("Utility function must be specified")
        if args.util not in ['split', 'binarize', 'validate_survival']:
            raise ValueError(f"Invalid utility function: {args.util}")

        if args.util == 'split':
            # Validate inputs
            if not args.clin:
                raise ValueError("Clinical data file must be provided")
            if not args.omics:
                raise ValueError("At least one omics file must be provided")
            if not os.path.isfile(args.clin):
                raise FileNotFoundError(f"Clinical file not found: {args.clin}")
            # Validate split ratio
            if not 0 < args.split < 1:
                raise ValueError(f"Split ratio must be between 0 and 1, got {args.split}")

        elif args.util == 'binarize':
            # Validate mutation data file
            if not args.mutations:
                raise ValueError("Mutation data file must be provided")
            if not os.path.isfile(args.mutations):
                raise FileNotFoundError(f"Mutation data file not found: {args.mutations}")
            # Validate gene and sample indices
            if args.gene_idx < 0 or args.sample_idx < 0:
                raise ValueError("Gene and sample indices must be non-negative integers")

        elif args.util == 'validate_survival':
            # Validate clinical data file
            if not args.clin:
                raise ValueError("Clinical data file must be provided")
            if not os.path.isfile(args.clin):
                raise FileNotFoundError(f"Clinical file not found: {args.clin}")
            # Validate survival event variable
            if not args.surv_event_var:
                raise ValueError("Survival event variable must be specified")

        # Create output directory if it doesn't exist
        if not os.path.exists(args.out):
            os.makedirs(args.out)

        if args.util == 'split':
            # Parse omics files
            omics_files = [f.strip() for f in args.omics.split(',') if f.strip()]
            if not omics_files:
                raise ValueError("At least one omics file must be provided")
            # Check omics files exist
            for f in omics_files:
                if not os.path.isfile(f):
                    raise FileNotFoundError(f"Omics file not found: {f}")
            data = make_data_dict(args.clin, omics_files)
            split_and_save_data(data, ratio=args.split, output_dir=args.out)

        elif args.util == 'binarize':
            mutations_df = read_data(args.mutations, index=False)
            if mutations_df.empty:
                raise ValueError("Mutation data file is empty")

            binarized_matrix = binarize_mutations(mutations_df, gene_idx=args.gene_idx, sample_idx=args.sample_idx)
            # Save binarized matrix
            output_file = os.path.join(args.out, 'binarized_mutations.tabular')
            binarized_matrix.to_csv(output_file, sep='\t')
            print(f"Binarized mutation matrix saved to {output_file}")

        elif args.util == 'validate_survival':
            clin_df = read_data(args.clin, index=False)
            if clin_df.empty:
                raise ValueError("Clinical data file is empty")

            # Validate survival event variable
            validate_survival(clin_df, args.surv_event_var)
            print(f"Survival event variable '{args.surv_event_var}' validated successfully.")

    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
