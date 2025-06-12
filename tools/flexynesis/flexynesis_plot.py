#!/usr/bin/env python
"""Generate plots using flexynesis
This script generates dimensionality reduction plots, Kaplan-Meier survival curves,
and Cox proportional hazards models from data processed by flexynesis."""

import argparse
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import torch
from flexynesis import build_cox_model, get_important_features, plot_dim_reduced, plot_hazard_ratios, plot_kaplan_meier_curves, plot_scatter


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


def load_omics(omics_path):
    """Load omics data from a file. First column should be features"""
    try:
        # Determine file extension
        file_ext = Path(omics_path).suffix.lower()

        if file_ext == '.csv':
            df = pd.read_csv(omics_path, index_col=0)
        elif file_ext in ['.tsv', '.txt', '.tab', '.tabular']:
            df = pd.read_csv(omics_path, sep='\t', index_col=0)
        else:
            raise ValueError(f"Unsupported file extension: {file_ext}")
        return df

    except Exception as e:
        raise ValueError(f"Error loading omics data from {omics_path}: {e}") from e


def load_model(model_path):
    """Load flexynesis model from pickle file"""
    try:
        with open(model_path, 'rb') as f:
            model = torch.load(f, weights_only=False)
        return model
    except Exception as e:
        raise ValueError(f"Error loading model from {model_path}: {e}") from e


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


def plot_label_concordance_heatmap(labels1, labels2, figsize=(12, 10)):
    """
    Plot a heatmap reflecting the concordance between two sets of labels using pandas crosstab.

    Parameters:
    - labels1: The first set of labels.
    - labels2: The second set of labels.
    """
    # Compute the cross-tabulation
    ct = pd.crosstab(pd.Series(labels1, name='Labels Set 1'), pd.Series(labels2, name='Labels Set 2'))
    # Normalize the cross-tabulation matrix column-wise
    ct_normalized = ct.div(ct.sum(axis=1), axis=0)

    # Plot the heatmap
    plt.figure(figsize=figsize)
    sns.heatmap(ct_normalized, annot=True, cmap='viridis', linewidths=.5)  # col_cluster=False)
    plt.title('Concordance between label groups')

    return plt.gcf()


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

    # Reset index and rename the index column to sample_id
    survival_data = survival_data.reset_index()
    if survival_data.columns[0] != 'sample_id':
        survival_data = survival_data.rename(columns={survival_data.columns[0]: 'sample_id'})

    # Convert survival event column to binary (0/1) based on event_value
    # Check if the event column exists
    if args.surv_event_var not in survival_data.columns:
        raise ValueError(f"Column '{args.surv_event_var}' not found in survival data")

    # Convert to string for comparison to handle mixed types
    survival_data[args.surv_event_var] = survival_data[args.surv_event_var].astype(str)
    event_value_str = str(args.event_value)

    # Create binary event column (1 if matches event_value, 0 otherwise)
    survival_data[f'{args.surv_event_var}_binary'] = (
        survival_data[args.surv_event_var] == event_value_str
    ).astype(int)

    # Filter for survival category and class_label == '1:DECEASED'
    label_data['class_label'] = label_data['class_label'].astype(str)

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
        events=df_deceased[f'{args.surv_event_var}_binary'],
        categorical_variable=group_labels
    )

    output_path_known = output_dir / f"{output_name_base}_km_risk_subtypes.{args.format}"
    print(f"Saving Kaplan-Meier plot to: {output_path_known.absolute()}")
    fig_known.save(output_path_known, dpi=args.dpi, bbox_inches='tight')

    print("Kaplan-Meier plot saved successfully!")


def generate_cox_plots(model, clinical_train, clinical_test, omics_train, omics_test, args, output_dir, output_name_base):
    """Generate Cox proportional hazards plots"""
    print("Generating Cox proportional hazards analysis...")

    # Parse clinical variables
    clinical_vars = [var.strip() for var in args.clinical_variables.split(',')]

    # Validate that survival variables are included
    required_vars = [args.surv_time_var, args.surv_event_var]
    for var in required_vars:
        if var not in clinical_vars:
            clinical_vars.append(var)

    print(f"Using clinical variables: {', '.join(clinical_vars)}")

    # filter datasets for clinical variables
    if all(var in clinical_train.columns and var in clinical_test.columns for var in clinical_vars):
        df_clin_train = clinical_train[clinical_vars]
        df_clin_test = clinical_test[clinical_vars]
        # Drop rows with NaN in clinical variables
        df_clin_train = df_clin_train.dropna(subset=clinical_vars)
        df_clin_test = df_clin_test.dropna(subset=clinical_vars)
    else:
        raise ValueError(f"Not all clinical variables found in datasets. Available in train dataset: {clinical_train.columns.tolist()}, Available in test dataset: {clinical_test.columns.tolist()}")

    # Combine
    df_clin = pd.concat([df_clin_train, df_clin_test], axis=0)

    # Get top survival markers
    print(f"Extracting top {args.top_features} important features for {args.surv_event_var}...")
    try:
        imp = get_important_features(model,
                                     var=args.surv_event_var,
                                     top=args.top_features
                                     )['name'].unique().tolist()
        print(f"Top features: {', '.join(imp)}")
    except Exception as e:
        raise ValueError(f"Error getting important features: {e}")

    # Extract feature data from omics datasets
    try:
        omics_test = omics_test.loc[omics_test.index.isin(imp)]
        omics_train = omics_train.loc[omics_train.index.isin(imp)]
        # Drop rows with NaN in omics datasets
        omics_test = omics_test.dropna(subset=omics_test.columns)
        omics_train = omics_train.dropna(subset=omics_train.columns)

        df_imp = pd.concat([omics_train, omics_test], axis=1)
        df_imp = df_imp.T  # Transpose to have samples as rows

        print(f"Feature data shape: {df_imp.shape}")
    except Exception as e:
        raise ValueError(f"Error extracting feature subset: {e}")

    # Combine markers with clinical variables
    df = pd.merge(df_imp, df_clin, left_index=True, right_index=True)
    print(f"Combined data shape: {df.shape}")

    # Remove samples without survival endpoints
    initial_samples = len(df)
    df = df[df[args.surv_event_var].notna()]
    final_samples = len(df)
    print(f"Removed {initial_samples - final_samples} samples without survival data")

    if df.empty:
        raise ValueError("No samples remain after filtering for survival data")

    # Convert survival event column to binary (0/1) based on event_value
    # Convert to string for comparison to handle mixed types
    df[args.surv_event_var] = df[args.surv_event_var].astype(str)
    event_value_str = str(args.event_value)

    df[f'{args.surv_event_var}'] = (
        df[args.surv_event_var] == event_value_str
    ).astype(int)

    # Build Cox model
    print(f"Building Cox model with time variable: {args.surv_time_var}, event variable: {args.surv_event_var}")
    try:
        coxm = build_cox_model(df,
                               duration_col=args.surv_time_var,
                               event_col=args.surv_event_var)
        print("Cox model built successfully")
    except Exception as e:
        raise ValueError(f"Error building Cox model: {e}")

    # Generate hazard ratios plot
    try:
        print("Generating hazard ratios plot...")
        fig = plot_hazard_ratios(coxm)

        output_path = output_dir / f"{output_name_base}_hazard_ratios.{args.format}"
        print(f"Saving hazard ratios plot to: {output_path.absolute()}")
        fig.save(output_path, dpi=args.dpi, bbox_inches='tight')

        print("Cox proportional hazards analysis completed successfully!")

    except Exception as e:
        raise ValueError(f"Error generating hazard ratios plot: {e}")


def generate_plot_scatter(labels, args, output_dir, output_name_base):
    """Generate scatter plot of known vs predicted labels"""
    print("Generating scatter plot of known vs predicted labels...")

    # filter labels for the target value
    if args.target_value:
        labels = labels[labels['variable'] == args.target_value]
    if labels.empty:
        raise ValueError(f"No data found for target value '{args.target_value}' in labels")

    true_values = pd.to_numeric(labels['known_label'], errors='coerce')
    predicted_values = pd.to_numeric(labels['predicted_label'], errors='coerce')
    if true_values.isna().all() or predicted_values.isna().all():
        raise ValueError("No valid numeric values found for known or predicted labels")

    print(f"Plotting scatter plot for target value: {args.target_value}")
    fig = plot_scatter(true_values, predicted_values)

    output_path = output_dir / f"{output_name_base}.{args.format}"
    print(f"Saving scatter plot to: {output_path.absolute()}")
    fig.save(output_path, dpi=args.dpi, bbox_inches='tight')

    print("Scatter plot generated successfully!")


def generate_label_concordance_heatmap(labels, args, output_dir, output_name_base):
    """Generate label concordance heatmap"""
    print("Generating label concordance heatmap...")

    # Filter labels for the target value
    if args.target_value:
        labels = labels[labels['variable'] == args.target_value]
    if labels.empty:
        raise ValueError(f"No data found for target value '{args.target_value}' in labels")

    true_values = labels['known_label'].tolist()
    predicted_values = labels['predicted_label'].tolist()

    print("Plotting label concordance heatmap...")
    fig = plot_label_concordance_heatmap(true_values, predicted_values)
    plt.close(fig)

    output_path = output_dir / f"{output_name_base}_heatmap.{args.format}"
    print(f"Saving heatmap to: {output_path.absolute()}")
    fig.savefig(output_path, dpi=args.dpi, bbox_inches='tight')

    print("Label concordance heatmap generated successfully!")


def main():
    """Main function to parse arguments and generate plots"""
    parser = argparse.ArgumentParser(description="Generate plots using flexynesis")

    # Required arguments
    parser.add_argument("--labels", type=str, required=False,
                        help="Path to labels file generated by flexynesis")

    # Plot type
    parser.add_argument("--plot_type", type=str, required=True,
                        choices=['dimred', 'kaplan_meier', 'cox', 'scatter', 'concordance_heatmap'],
                        help="Type of plot to generate: 'dimred' for dimensionality reduction, 'kaplan_meier' for survival analysis, 'cox' for Cox proportional hazards")

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

    # Arguments for Cox analysis
    parser.add_argument("--model", type=str,
                        help="Path to trained flexynesis model (pickle file). Required for cox plots.")
    parser.add_argument("--clinical_train", type=str,
                        help="Path to training dataset (pickle file). Required for cox plots.")
    parser.add_argument("--clinical_test", type=str,
                        help="Path to test dataset (pickle file). Required for cox plots.")
    parser.add_argument("--omics_train", type=str, default=None,
                        help="Path to training omics dataset. Optional for cox plots.")
    parser.add_argument("--omics_test", type=str, default=None,
                        help="Path to test omics dataset. Optional for cox plots.")
    parser.add_argument("--clinical_variables", type=str,
                        help="Comma-separated list of clinical variables to include in Cox model (e.g., 'AGE,SEX,HISTOLOGICAL_DIAGNOSIS,STUDY')")
    parser.add_argument("--top_features", type=int, default=20,
                        help="Number of top important features to include in Cox model. Default is 5")

    # Arguments for scatter plot and heatmap
    parser.add_argument("--target_value", type=str, default=None,
                        help="Target value for scatter plot.")

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
        # validate plot type
        if not args.plot_type:
            raise ValueError("Please specify a plot type using --plot_type")
        if args.plot_type not in ['dimred', 'kaplan_meier', 'cox', 'scatter', 'concordance_heatmap']:
            raise ValueError(f"Invalid plot type: {args.plot_type}. Must be one of: 'dimred', 'kaplan_meier', 'cox', 'scatter', 'concordance_heatmap'")

        # Validate plot type requirements
        if args.plot_type in ['dimred']:
            if not args.embeddings:
                raise ValueError("--embeddings is required when plot_type is 'dimred'")
            if not os.path.isfile(args.embeddings):
                raise FileNotFoundError(f"embeddings file not found: {args.embeddings}")
            if not args.labels:
                raise ValueError("--labels is required for dimensionality reduction plots")
            if not args.method:
                raise ValueError("--method is required for dimensionality reduction plots")
            if not args.target_variables:
                raise ValueError("--target_variables is required for dimensionality reduction plots")

        if args.plot_type in ['kaplan_meier']:
            if not args.survival_data:
                raise ValueError("--survival_data is required when plot_type is 'kaplan_meier'")
            if not os.path.isfile(args.survival_data):
                raise FileNotFoundError(f"Survival data file not found: {args.survival_data}")
            if not args.labels:
                raise ValueError("--labels is required for dimensionality reduction plots")
            if not args.method:
                raise ValueError("--method is required for dimensionality reduction plots")
            if not args.surv_time_var:
                raise ValueError("--surv_time_var is required for Kaplan-Meier plots")
            if not args.surv_event_var:
                raise ValueError("--surv_event_var is required for Kaplan-Meier plots")
            if not args.event_value:
                raise ValueError("--event_value is required for Kaplan-Meier plots")

        if args.plot_type in ['cox']:
            if not args.model:
                raise ValueError("--model is required when plot_type is 'cox'")
            if not os.path.isfile(args.model):
                raise FileNotFoundError(f"Model file not found: {args.model}")
            if not args.clinical_train:
                raise ValueError("--clinical_train is required when plot_type is 'cox'")
            if not os.path.isfile(args.clinical_train):
                raise FileNotFoundError(f"Training dataset file not found: {args.clinical_train}")
            if not args.clinical_test:
                raise ValueError("--clinical_test is required when plot_type is 'cox'")
            if not os.path.isfile(args.clinical_test):
                raise FileNotFoundError(f"Test dataset file not found: {args.clinical_test}")
            if not args.omics_train:
                raise ValueError("--omics_train is required when plot_type is 'cox'")
            if not os.path.isfile(args.omics_train):
                raise FileNotFoundError(f"Training omics dataset file not found: {args.omics_train}")
            if not args.omics_test:
                raise ValueError("--omics_test is required when plot_type is 'cox'")
            if not os.path.isfile(args.omics_test):
                raise FileNotFoundError(f"Test omics dataset file not found: {args.omics_test}")
            if not args.surv_time_var:
                raise ValueError("--surv_time_var is required for Cox plots")
            if not args.surv_event_var:
                raise ValueError("--surv_event_var is required for Cox plots")
            if not args.clinical_variables:
                raise ValueError("--clinical_variables is required for Cox plots")
            if not isinstance(args.top_features, int) or args.top_features <= 0:
                raise ValueError("--top_features must be a positive integer")
            if not args.event_value:
                raise ValueError("--event_value is required for Kaplan-Meier plots")

        if args.plot_type in ['scatter']:
            if not args.labels:
                raise ValueError("--labels is required for scatter plots")
            if not args.target_value:
                raise ValueError("--target_value is required for scatter plots")
            if not os.path.isfile(args.labels):
                raise FileNotFoundError(f"Labels file not found: {args.labels}")

        if args.plot_type in ['concordance_heatmap']:
            if not args.labels:
                raise ValueError("--labels is required for concordance heatmap")
            if not args.target_value:
                raise ValueError("--target_value is required for concordance heatmap")
            if not os.path.isfile(args.labels):
                raise FileNotFoundError(f"Labels file not found: {args.labels}")

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
            elif args.plot_type == 'cox':
                model_name = Path(args.model).stem
                output_name_base = f"{model_name}_cox"
            elif args.plot_type == 'scatter':
                labels_name = Path(args.labels).stem
                output_name_base = f"{labels_name}_scatter"
            elif args.plot_type == 'concordance_heatmap':
                labels_name = Path(args.labels).stem
                output_name_base = f"{labels_name}_concordance"

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

        elif args.plot_type in ['cox']:
            # Load model and datasets
            print(f"Loading model from: {args.model}")
            model = load_model(args.model)
            print(f"Loading training dataset from: {args.clinical_train}")
            clinical_train = load_omics(args.clinical_train)
            print(f"Loading test dataset from: {args.clinical_test}")
            clinical_test = load_omics(args.clinical_test)
            print(f"Loading training omics dataset from: {args.omics_train}")
            omics_train = load_omics(args.omics_train)
            print(f"Loading test omics dataset from: {args.omics_test}")
            omics_test = load_omics(args.omics_test)

            generate_cox_plots(model, clinical_train, clinical_test, omics_test, omics_train, args, output_dir, output_name_base)

        elif args.plot_type in ['scatter']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)

            generate_plot_scatter(label_data, args, output_dir, output_name_base)

        elif args.plot_type in ['concordance_heatmap']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)

            generate_label_concordance_heatmap(label_data, args, output_dir, output_name_base)

        print("All plots generated successfully!")

    except (FileNotFoundError, ValueError, pd.errors.ParserError) as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
