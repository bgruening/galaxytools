#!/usr/bin/env python
"""Generate plots using flexynesis
This script generates dimensionality reduction plots and Kaplan-Meier survival curves
from data processed by flexynesis."""

import argparse
import os
from pathlib import Path

import pandas as pd
import numpy as np
from flexynesis import plot_dim_reduced, plot_kaplan_meier_curves


def load_embeddings(embeddings_path):
    """Load embeddings from a file"""
    try:
        # Determine file extension
        file_ext = Path(embeddings_path).suffix.lower()

        if file_ext == '.csv':
            df = pd.read_csv(embeddings_path, index_col=0)
        elif file_ext in ['.tsv', '.txt', '.tab', '.tabular']:
            df = pd.read_csv(embeddings_path, sep='\t', index_col=0)
        else:
            raise ValueError(f"Unsupported file extension: {file_ext}")

        return df, df.index.tolist()

    except Exception as e:
        raise ValueError(f"Error loading embeddings from {embeddings_path}: {e}") from e


def load_labels(labels_input, plot_type=None):
    """Load predicted labels from flexynesis"""
    try:
        # Determine file extension
        file_ext = Path(labels_input).suffix.lower()

        if file_ext == '.csv':
            df = pd.read_csv(labels_input)
        elif file_ext in ['.tsv', '.txt', '.tab', '.tabular']:
            df = pd.read_csv(labels_input, sep='\t')

        # Check if this is the specific format with sample_id, known_label, predicted_label
        required_cols = ['sample_id', 'variable', 'class_label', 'probability', 'known_label', 'predicted_label']
        if all(col in df.columns for col in required_cols):
            return df
        else:
            raise ValueError(f"Labels file {labels_input} does not contain required columns: {required_cols}")

    except Exception as e:
        raise ValueError(f"Error loading labels from {labels_input}: {e}") from e


def load_survival_data(survival_path):
    """Load survival data from a file. First column should be sample_id"""
    try:
        # Determine file extension
        file_ext = Path(survival_path).suffix.lower()

        if file_ext == '.csv':
            df = pd.read_csv(survival_path, index_col=0)
        elif file_ext in ['.tsv', '.txt', '.tab', '.tabular']:
            df = pd.read_csv(survival_path, sep='\t', index_col=0)
        else:
            raise ValueError(f"Unsupported file extension: {file_ext}")
        return df

    except Exception as e:
        raise ValueError(f"Error loading survival data from {survival_path}: {e}") from e


def match_samples_to_embeddings(sample_names, label_data):
    """Filter label data to match sample names in the embeddings"""
    df_matched = label_data[label_data['sample_id'].isin(sample_names)]
    return df_matched


def detect_color_type(labels_series):
    """Auto-detect whether target variables should be treated as categorical or numerical"""
    # Remove NaN
    clean_labels = labels_series.dropna()

    if clean_labels.empty:
        return 'categorical'  # default output if no labels

    # Check if all values can be converted to numbers
    try:
        numeric_labels = pd.to_numeric(clean_labels, errors='coerce')

        # If conversion failed -> categorical
        if numeric_labels.isna().any():
            return 'categorical'

        # Check number of unique values
        unique_count = len(clean_labels.unique())
        total_count = len(clean_labels)

        # If few unique values relative to total -> categorical
        # Threshold: if unique values < 10 OR unique/total < 0.1
        if unique_count < 10 or (unique_count / total_count) < 0.1:
            return 'categorical'
        else:
            return 'numerical'

    except Exception:
        return 'categorical'


def generate_dimred_plots(embeddings, matched_labels, args, output_dir, output_name_base):
    """Generate dimensionality reduction plots"""

    # Parse target variables
    target_vars = [var.strip() for var in args.target_variables.split(',')]

    print(f"Generating {args.method.upper()} plots for {len(target_vars)} target variable(s): {', '.join(target_vars)}")

    # Check variables
    available_vars = matched_labels['variable'].unique()
    missing_vars = [var for var in target_vars if var not in available_vars]

    if missing_vars:
        print(f"Warning: The following target variables were not found in the data: {', '.join(missing_vars)}")
        print(f"Available variables: {', '.join(available_vars)}")

    # Filter to only process available variables
    valid_vars = [var for var in target_vars if var in available_vars]

    if not valid_vars:
        raise ValueError(f"None of the specified target variables were found in the data. Available: {', '.join(available_vars)}")

    # Generate plots for each valid target variable
    for var in valid_vars:
        print(f"\nPlotting variable: {var}")

        # Filter matched labels for current variable
        var_labels = matched_labels[matched_labels['variable'] == var].copy()
        var_labels = var_labels.drop_duplicates(subset='sample_id')

        if var_labels.empty:
            print(f"Warning: No data found for variable '{var}', skipping...")
            continue

        # Auto-detect color type
        known_color_type = detect_color_type(var_labels['known_label'])
        predicted_color_type = detect_color_type(var_labels['predicted_label'])

        print(f"  Auto-detected color types - Known: {known_color_type}, Predicted: {predicted_color_type}")

        try:
            # Plot 1: Known labels
            print(f"  Creating known labels plot for {var}...")
            fig_known = plot_dim_reduced(
                matrix=embeddings,
                labels=var_labels['known_label'],
                method=args.method,
                color_type=known_color_type
            )

            output_path_known = output_dir / f"{output_name_base}_{var}_known.{args.format}"
            print(f"  Saving known labels plot to: {output_path_known.name}")
            fig_known.save(output_path_known, dpi=args.dpi, bbox_inches='tight')

            # Plot 2: Predicted labels
            print(f"  Creating predicted labels plot for {var}...")
            fig_predicted = plot_dim_reduced(
                matrix=embeddings,
                labels=var_labels['predicted_label'],
                method=args.method,
                color_type=predicted_color_type
            )

            output_path_predicted = output_dir / f"{output_name_base}_{var}_predicted.{args.format}"
            print(f"  Saving predicted labels plot to: {output_path_predicted.name}")
            fig_predicted.save(output_path_predicted, dpi=args.dpi, bbox_inches='tight')

            print(f"  ✓ Successfully created plots for variable '{var}'")

        except Exception as e:
            print(f"  ✗ Error creating plots for variable '{var}': {e}")
            continue

    print(f"\nDimensionality reduction plots completed for {len(valid_vars)} variable(s)!")


def generate_km_plots(survival_data, label_data, args, output_dir, output_name_base):
    """Generate Kaplan-Meier plots"""
    print("Generating Kaplan-Meier curves of risk subtypes...")

    survival_data = survival_data.reset_index().rename(columns={'index': 'sample_id'})

    # Filter for survival category and class_label == '1:DECEASED'
    label_data['class_label'] = label_data['class_label'].astype(str)
    # Convert args.event_value to string for consistent comparison
    event_value_str = str(args.event_value)

    label_data = label_data[(label_data['variable'] == args.surv_event_var) & (label_data['class_label'] == event_value_str)]

    # check survival data
    for col in [args.surv_time_var, args.surv_event_var]:
        if col not in survival_data.columns:
            raise ValueError(f"Column '{col}' not found in survival data")

    # Merge survival data with labels
    df_deceased = pd.merge(survival_data, label_data, on='sample_id', how='inner')

    if df_deceased.empty:
        raise ValueError("No matching samples found after merging survival and label data.")

    # Get risk scores
    risk_scores = df_deceased['probability'].values

    # Compute groups (e.g., median split)
    quantiles = np.quantile(risk_scores, [0.5])
    groups = np.digitize(risk_scores, quantiles)
    group_labels = ['low_risk' if g == 0 else 'high_risk' for g in groups]

    fig_known = plot_kaplan_meier_curves(
        durations=df_deceased[args.surv_time_var],
        events=df_deceased[args.surv_event_var],
        categorical_variable=group_labels
    )

    output_path_known = output_dir / f"{output_name_base}_km_risk_subtypes.{args.format}"
    print(f"Saving Kaplan-Meier plot to: {output_path_known.absolute()}")
    fig_known.save(output_path_known, dpi=args.dpi, bbox_inches='tight')

    print("Kaplan-Meier plot saved successfully!")


def main():
    """Main function to parse arguments and generate plots"""
    parser = argparse.ArgumentParser(description="Generate plots using flexynesis")

    # Required arguments
    parser.add_argument("--labels", type=str, required=True,
                        help="Path to labels file generated by flexynesis")

    # Plot type
    parser.add_argument("--plot_type", type=str, required=True,
                        choices=['dimred', 'kaplan_meier'],
                        help="Type of plot to generate: 'dimred' for dimensionality reduction, 'kaplan_meier' for survival analysis")

    # Arguments for dimensionality reduction
    parser.add_argument("--embeddings", type=str,
                        help="Path to input data embeddings file (CSV or tabular format). Required for dimred plots.")
    parser.add_argument("--method", type=str, default='pca', choices=['pca', 'umap'],
                        help="Transformation method ('pca' or 'umap'). Default is 'pca'. Used for dimred plots.")
    parser.add_argument("--target_variables", type=str, required=False,
                        help="Comma-separated list of target variables to plot.")

    # Arguments for Kaplan-Meier
    parser.add_argument("--survival_data", type=str,
                        help="Path to survival data file with columns: duration and event. Required for kaplan_meier plots.")
    parser.add_argument("--surv_time_var", type=str, required=False,
                        help="Column name for survival time")
    parser.add_argument("--surv_event_var", type=str, required=False,
                        help="Column name for survival event")
    parser.add_argument("--event_value", type=str, required=False,
                        help="Value in event column that represents an event (e.g., 'DECEASED')")

    # Common arguments
    parser.add_argument("--output_dir", type=str, default='output',
                        help="Output directory. Default is 'output'")
    parser.add_argument("--output_name", type=str, default=None,
                        help="Output filename base")
    parser.add_argument("--format", type=str, default='jpg', choices=['png', 'pdf', 'svg', 'jpg'],
                        help="Output format for the plot. Default is 'jpg'")
    parser.add_argument("--dpi", type=int, default=300,
                        help="DPI for the output image. Default is 300")

    args = parser.parse_args()

    try:
        # Validate plot type requirements
        if args.plot_type in ['dimred']:
            if not args.embeddings:
                raise ValueError("--embeddings is required when plot_type is 'dimred'")
            if not os.path.isfile(args.embeddings):
                raise FileNotFoundError(f"embeddings file not found: {args.embeddings}")
            if not args.target_variables:
                raise ValueError("--target_variables is required for dimensionality reduction plots")

        if args.plot_type in ['kaplan_meier']:
            if not args.survival_data:
                raise ValueError("--survival_data is required when plot_type is 'kaplan_meier'")
            if not os.path.isfile(args.survival_data):
                raise FileNotFoundError(f"Survival data file not found: {args.survival_data}")
            if not args.surv_time_var:
                raise ValueError("--surv_time_var is required for Kaplan-Meier plots")
            if not args.surv_event_var:
                raise ValueError("--surv_event_var is required for Kaplan-Meier plots")
            if not args.event_value:
                raise ValueError("--event_value is required for Kaplan-Meier plots")

        # Validate other arguments
        if args.method not in ['pca', 'umap']:
            raise ValueError("Method must be 'pca' or 'umap'")

        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        print(f"Output directory: {output_dir.absolute()}")

        # Generate output filename base
        if args.output_name:
            output_name_base = args.output_name
        else:
            if args.plot_type == 'dimred':
                embeddings_name = Path(args.embeddings).stem
                output_name_base = f"{embeddings_name}_{args.method}"
            elif args.plot_type == 'kaplan_meier':
                survival_name = Path(args.survival_data).stem
                output_name_base = f"{survival_name}_km"

        # Generate plots based on type
        if args.plot_type in ['dimred']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels, plot_type='dimred')
            # Load embeddings data
            print(f"Loading embeddings from: {args.embeddings}")
            embeddings, sample_names = load_embeddings(args.embeddings)
            print(f"embeddings shape: {embeddings.shape}")

            # Match samples to embeddings
            matched_labels = match_samples_to_embeddings(sample_names, label_data)
            print(f"Successfully matched {len(matched_labels)} samples for dimensionality reduction")

            generate_dimred_plots(embeddings, matched_labels, args, output_dir, output_name_base)

        elif args.plot_type in ['kaplan_meier']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)
            # Load survival data
            print(f"Loading survival data from: {args.survival_data}")
            survival_data = load_survival_data(args.survival_data)
            print(f"Survival data shape: {survival_data.shape}")

            generate_km_plots(survival_data, label_data, args, output_dir, output_name_base)

        print("All plots generated successfully!")

    except (FileNotFoundError, ValueError, pd.errors.ParserError) as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
