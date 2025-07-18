#!/usr/bin/env python
"""Generate plots using flexynesis
This script generates dimensionality reduction plots, Kaplan-Meier survival curves,
and Cox proportional hazards models from data processed by flexynesis."""

import argparse
import os
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from flexynesis import (
    build_cox_model,
    plot_dim_reduced,
    plot_hazard_ratios,
    plot_kaplan_meier_curves,
    plot_pr_curves,
    plot_roc_curves,
    plot_scatter
)
from scipy.stats import kruskal, mannwhitneyu


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


def load_labels(labels_input):
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
        print(f"available columns: {df.columns.tolist()}")
        if all(col in df.columns for col in required_cols):
            print("Detected flexynesis labels format")
        else:
            print("Labels are not in flexynesis format (Custom labels)")

        return df

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


def match_samples_to_embeddings(sample_names, label_data):
    """Filter label data to match sample names in the embeddings"""
    # Create a DataFrame from sample_names to preserve order
    sample_df = pd.DataFrame({'sample_names': sample_names})

    # left_join
    first_column = label_data.columns[0]
    df_matched = sample_df.merge(label_data, left_on='sample_names', right_on=first_column, how='left')

    # remove sample_names to keep the initial structure
    df_matched = df_matched.drop('sample_names', axis=1)
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


def plot_boxplot(categorical_x, numerical_y, title_x='Categories', title_y='Values', figsize=(10, 6), jittersize=4):
    """
    Create a boxplot with to visualize the distribution of predicted probabilities across different categories.
    the x axis represents the true labels, and the y axis represents the predicted probabilities for specific categories.
    """
    df = pd.DataFrame({title_x: categorical_x, title_y: numerical_y})

    # Compute p-value
    groups = df[title_x].unique()
    if len(groups) == 2:
        group1 = df[df[title_x] == groups[0]][title_y]
        group2 = df[df[title_x] == groups[1]][title_y]
        stat, p = mannwhitneyu(group1, group2, alternative='two-sided')
        test_name = "Mann-Whitney U"
    else:
        group_data = [df[df[title_x] == group][title_y] for group in groups]
        stat, p = kruskal(*group_data)
        test_name = "Kruskal-Wallis"

    # Create a boxplot with jittered points
    plt.figure(figsize=figsize)
    sns.boxplot(x=title_x, y=title_y, hue=title_x, data=df, palette='Set2', legend=False, fill=False)
    sns.stripplot(x=title_x, y=title_y, data=df, color='black', size=jittersize, jitter=True, dodge=True, alpha=0.4)

    # Labels and p-value annotation
    plt.xlabel(title_x)
    plt.ylabel(title_y)
    plt.text(
        x=-0.4,
        y=plt.ylim()[1],
        s=f'{test_name} p = {p:.3e}',
        verticalalignment='top',
        horizontalalignment='left',
        fontsize=12,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='gray')
    )

    plt.tight_layout()
    return plt.gcf()


def generate_dimred_plots(embeddings, matched_labels, args, output_dir, output_name_base):
    """Generate dimensionality reduction plots"""

    # Check if this is the specific format with sample_id, known_label, predicted_label
    required_cols = ['sample_id', 'variable', 'class_label', 'probability', 'known_label', 'predicted_label']
    if all(col in matched_labels.columns for col in required_cols):
        print("Detected flexynesis labels format")
        # Parse target values from comma-separated string
        if args.target_value:
            target_values = [val.strip() for val in args.target_value.split(',')]
        else:
            # If no target values specified, use all unique variables
            target_values = matched_labels['variable'].unique().tolist()

        print(f"Generating {args.method.upper()} plots for {len(target_values)} target variable(s): {', '.join(target_values)}")

        # Check variables
        available_vars = matched_labels['variable'].unique()
        missing_vars = [var for var in target_values if var not in available_vars]

        if missing_vars:
            print(f"Warning: The following target variables were not found in the data: {', '.join(missing_vars)}")
            print(f"Available variables: {', '.join(available_vars)}")

        # Filter to only process available variables
        valid_vars = [var for var in target_values if var in available_vars]

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

    else:
        print("Labels are not in flexynesis format (Custom labels)")
        print(f"Generating {args.method.upper()} plots")

        # check if the color argument is provided
        if args.color is None:
            raise ValueError("No color argument provided. Please specify a color variable.")

        # check if the color variable exists in matched_labels
        if args.color not in matched_labels.columns:
            raise ValueError(f"Color variable '{args.color}' not found in matched labels. Available columns: {matched_labels.columns.tolist()}")

        # Auto-detect color type
        color_type = detect_color_type(matched_labels[args.color])

        print(f"  Auto-detected color types - Known: {color_type}")

        # Plot: Known labels
        print(f"  Creating known labels plot for {args.color}...")
        fig_known = plot_dim_reduced(
            matrix=embeddings,
            labels=matched_labels[args.color],
            method=args.method,
            color_type=color_type
        )

        output_path_known = output_dir / f"{output_name_base}_{args.color}_known.{args.format}"
        print(f"  Saving known labels plot to: {output_path_known.name}")
        fig_known.save(output_path_known, dpi=args.dpi, bbox_inches='tight')

        print(f"  ✓ Successfully created plots for variable '{args.color}'")


def generate_km_plots(survival_data, label_data, args, output_dir, output_name_base):
    """Generate Kaplan-Meier plots"""
    print("Generating Kaplan-Meier curves of risk subtypes...")

    # Reset index and rename the index column to sample_id
    survival_data = survival_data.reset_index()
    if survival_data.columns[0] != 'sample_id':
        survival_data = survival_data.rename(columns={survival_data.columns[0]: 'sample_id'})

    # Check if the event column exists
    if args.surv_event_var not in survival_data.columns:
        raise ValueError(f"Column '{args.surv_event_var}' not found in survival data")

    label_data = label_data[(label_data['variable'] == args.surv_event_var)]

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


def generate_cox_plots(important_features, clinical_train, clinical_test, omics_train, omics_test, args, output_dir, output_name_base):
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
        print(f"Loading {args.top_features} important features from: {args.important_features}")
        imp_features = load_labels(args.important_features)
        imp_features = imp_features[imp_features['target_variable'] == args.surv_event_var]
        if imp_features.empty:
            raise ValueError(f"No important features found for target variable '{args.surv_event_var}' in {args.important_features}")
        imp_features = imp_features.sort_values(by='importance', ascending=False)

        if len(imp_features) < args.top_features:
            raise ValueError(f"Requested top {args.top_features} features, but only {len(imp_features)} available in {args.important_features}")

        imp = imp_features['name'].unique().tolist()[0:args.top_features]

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

    # Build Cox model
    print(f"Building Cox model with time variable: {args.surv_time_var}, event variable: {args.surv_event_var}")
    try:
        coxm = build_cox_model(df,
                               duration_col=args.surv_time_var,
                               event_col=args.surv_event_var,
                               crossval=args.crossval,
                               n_splits=args.n_splits,
                               random_state=args.random_state)
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
    print("Generating scatter plots of known vs predicted labels...")

    # Parse target values from comma-separated string
    if args.target_value:
        target_values = [val.strip() for val in args.target_value.split(',')]
    else:
        # If no target values specified, use all unique variables
        target_values = labels['variable'].unique().tolist()

    print(f"Processing target values: {target_values}")

    successful_plots = 0
    skipped_plots = 0

    for target_value in target_values:
        print(f"\nProcessing target value: '{target_value}'")

        # Filter labels for the current target value
        target_labels = labels[labels['variable'] == target_value]

        if target_labels.empty:
            print(f"  Warning: No data found for target value '{target_value}' - skipping")
            skipped_plots += 1
            continue

        # Check if labels are numeric and convert
        true_values = pd.to_numeric(target_labels['known_label'], errors='coerce')
        predicted_values = pd.to_numeric(target_labels['predicted_label'], errors='coerce')

        if true_values.isna().all() or predicted_values.isna().all():
            print(f"No valid numeric values found for known or predicted labels in '{target_value}'")
            skipped_plots += 1
            continue

        try:
            print(f"  Generating scatter plot for '{target_value}'...")
            fig = plot_scatter(true_values, predicted_values)

            # Create output filename with target value
            safe_target_name = target_value.replace('/', '_').replace('\\', '_').replace(' ', '_')
            if len(target_values) > 1:
                output_filename = f"{output_name_base}_{safe_target_name}.{args.format}"
            else:
                output_filename = f"{output_name_base}.{args.format}"

            output_path = output_dir / output_filename
            print(f"  Saving scatter plot to: {output_path.absolute()}")
            fig.save(output_path, dpi=args.dpi, bbox_inches='tight')

            successful_plots += 1
            print(f"  Scatter plot for '{target_value}' generated successfully!")

        except Exception as e:
            print(f"  Error generating plot for '{target_value}': {str(e)}")
            skipped_plots += 1

    # Summary
    print("  Summary:")
    print(f"  Successfully generated: {successful_plots} plots")
    print(f"  Skipped: {skipped_plots} plots")

    if successful_plots == 0:
        raise ValueError("No scatter plots could be generated. Check your data and target values.")

    print("Scatter plot generation completed!")


def generate_label_concordance_heatmap(labels, args, output_dir, output_name_base):
    """Generate label concordance heatmap"""
    print("Generating label concordance heatmaps...")

    # Parse target values from comma-separated string
    if args.target_value:
        target_values = [val.strip() for val in args.target_value.split(',')]
    else:
        # If no target values specified, use all unique variables
        target_values = labels['variable'].unique().tolist()

    print(f"Processing target values: {target_values}")

    for target_value in target_values:
        print(f"\nProcessing target value: '{target_value}'")

        # Filter labels for the current target value
        target_labels = labels[labels['variable'] == target_value]

        if target_labels.empty:
            print(f"  Warning: No data found for target value '{target_value}' - skipping")
            continue

        true_values = target_labels['known_label'].tolist()
        predicted_values = target_labels['predicted_label'].tolist()

        try:
            print(f"  Generating heatmap for '{target_value}'...")
            fig = plot_label_concordance_heatmap(true_values, predicted_values)
            plt.close(fig)

            # Create output filename with target value
            safe_target_name = target_value.replace('/', '_').replace('\\', '_').replace(' ', '_')
            if len(target_values) > 1:
                output_filename = f"{output_name_base}_{safe_target_name}.{args.format}"
            else:
                output_filename = f"{output_name_base}.{args.format}"

            output_path = output_dir / output_filename
            print(f"  Saving heatmap to: {output_path.absolute()}")
            fig.savefig(output_path, dpi=args.dpi, bbox_inches='tight')

        except Exception as e:
            print(f"  Error generating heatmap for '{target_value}': {str(e)}")
            continue

    print("Label concordance heatmap generated successfully!")


def generate_pr_curves(labels, args, output_dir, output_name_base):
    """Generate precision-recall curves"""
    print("Generating precision-recall curves...")

    # Parse target values from comma-separated string
    if args.target_value:
        target_values = [val.strip() for val in args.target_value.split(',')]
    else:
        # If no target values specified, use all unique variables
        target_values = labels['variable'].unique().tolist()

    print(f"Processing target values: {target_values}")

    for target_value in target_values:
        print(f"\nProcessing target value: '{target_value}'")

        # Filter labels for the current target value
        target_labels = labels[labels['variable'] == target_value]

        # Check if this is a regression problem (no class probabilities)
        prob_columns = target_labels['class_label'].unique()
        non_na_probs = target_labels['probability'].notna().sum()

        print(f"  Class labels found: {list(prob_columns)}")
        print(f"  Non-NaN probabilities: {non_na_probs}/{len(target_labels)}")

        # If most probabilities are NaN, this is likely a regression problem
        if non_na_probs < len(target_labels) * 0.1:  # Less than 10% valid probabilities
            print("  Detected regression problem - precision-recall curves not applicable")
            print(f"  Skipping '{target_value}' (use regression evaluation metrics instead)")
            continue

        # Debug: Check data quality
        total_rows = len(target_labels)
        missing_labels = target_labels['known_label'].isna().sum()
        missing_probs = target_labels['probability'].isna().sum()
        unique_samples = target_labels['sample_id'].nunique()
        unique_classes = target_labels['class_label'].nunique()

        print(f"  Data summary: {total_rows} total rows, {unique_samples} unique samples, {unique_classes} unique classes")
        print(f"  Missing data: {missing_labels} missing known_label, {missing_probs} missing probability")

        if missing_labels > 0:
            print(f"  Warning: Found {missing_labels} missing known_label values")
            missing_samples = target_labels[target_labels['known_label'].isna()]['sample_id'].unique()[:5]
            print(f"  Sample IDs with missing known_label: {list(missing_samples)}")

            # Remove rows with missing known_label
            target_labels = target_labels.dropna(subset=['known_label'])
            if target_labels.empty:
                print(f"  Error: No valid known_label data remaining for '{target_value}' - skipping")
                continue

        # 1. Pivot to wide format
        prob_df = target_labels.pivot(index='sample_id', columns='class_label', values='probability')

        print(f"  After pivot: {prob_df.shape[0]} samples x {prob_df.shape[1]} classes")
        print(f"  Class columns: {list(prob_df.columns)}")

        # Check for NaN values in probability data
        nan_counts = prob_df.isna().sum()
        if nan_counts.any():
            print(f"  NaN counts per class: {dict(nan_counts)}")
            print(f"  Samples with any NaN: {prob_df.isna().any(axis=1).sum()}/{len(prob_df)}")

            # Drop only rows where ALL probabilities are NaN
            all_nan_rows = prob_df.isna().all(axis=1)
            if all_nan_rows.any():
                print(f"  Dropping {all_nan_rows.sum()} samples with all NaN probabilities")
                prob_df = prob_df[~all_nan_rows]

            remaining_nans = prob_df.isna().sum().sum()
            if remaining_nans > 0:
                print(f"  Warning: {remaining_nans} individual NaN values remain - filling with 0")
                prob_df = prob_df.fillna(0)

            if prob_df.empty:
                print(f"  Error: No valid probability data remaining for '{target_value}' - skipping")
                continue

        # 2. Get true labels
        true_labels_df = target_labels.drop_duplicates('sample_id')[['sample_id', 'known_label']].set_index('sample_id')

        # 3. Align indices - only keep samples that exist in both datasets
        common_indices = prob_df.index.intersection(true_labels_df.index)
        if len(common_indices) == 0:
            print(f"  Error: No common sample_ids between probability and true label data for '{target_value}' - skipping")
            continue

        print(f"  Found {len(common_indices)} samples with both probability and true label data")

        # Filter both datasets to common indices
        prob_df_aligned = prob_df.loc[common_indices]
        y_true = true_labels_df.loc[common_indices]['known_label']

        # 4. Final check for NaN values
        if y_true.isna().any():
            print(f"  Error: True labels still contain NaN after alignment for '{target_value}' - skipping")
            continue

        if prob_df_aligned.isna().any().any():
            print(f"  Error: Probability data still contains NaN after alignment for '{target_value}' - skipping")
            continue

        # 5. Convert categorical labels to integer labels
        # Create a mapping from class names to integers
        class_names = list(prob_df_aligned.columns)
        class_to_int = {class_name: i for i, class_name in enumerate(class_names)}

        print(f"  Class mapping: {class_to_int}")

        # Convert true labels to integers
        y_true_np = y_true.map(class_to_int).to_numpy()
        y_probs_np = prob_df_aligned.to_numpy()

        print(f"  Data shape: y_true={y_true_np.shape}, y_probs={y_probs_np.shape}")
        print(f"  Unique true labels (integers): {set(y_true_np)}")
        print(f"  Class labels (columns): {class_names}")
        print(f"  Label distribution: {dict(zip(*np.unique(y_true_np, return_counts=True)))}")

        # Check for any unmapped labels (will be NaN)
        if pd.isna(y_true_np).any():
            print("  Error: Some true labels could not be mapped to class columns")
            unmapped_labels = set(y_true[y_true.map(class_to_int).isna()])
            print(f"  Unmapped labels: {unmapped_labels}")
            print(f"  Available classes: {class_names}")
            continue

        try:
            print(f"  Generating precision-recall curve for '{target_value}'...")
            fig = plot_pr_curves(y_true_np, y_probs_np)

            # Create output filename with target value
            safe_target_name = target_value.replace('/', '_').replace('\\', '_').replace(' ', '_')
            if len(target_values) > 1:
                output_filename = f"{output_name_base}_{safe_target_name}.{args.format}"
            else:
                output_filename = f"{output_name_base}.{args.format}"

            output_path = output_dir / output_filename
            print(f"  Saving precision-recall curve to: {output_path.absolute()}")
            fig.save(output_path, dpi=args.dpi, bbox_inches='tight')

        except Exception as e:
            print(f"  Error generating precision-recall curve for '{target_value}': {str(e)}")
            print(f"  Debug info - y_true type: {type(y_true_np)}, contains NaN: {pd.isna(y_true_np).any()}")
            print(f"  Debug info - y_probs type: {type(y_probs_np)}, contains NaN: {pd.isna(y_probs_np).any()}")
            continue

    print("Precision-recall curves generated successfully!")


def generate_roc_curves(labels, args, output_dir, output_name_base):
    """Generate ROC curves"""
    print("Generating ROC curves...")

    # Parse target values from comma-separated string
    if args.target_value:
        target_values = [val.strip() for val in args.target_value.split(',')]
    else:
        # If no target values specified, use all unique variables
        target_values = labels['variable'].unique().tolist()

    print(f"Processing target values: {target_values}")

    for target_value in target_values:
        print(f"\nProcessing target value: '{target_value}'")

        # Filter labels for the current target value
        target_labels = labels[labels['variable'] == target_value]

        # Check if this is a regression problem (no class probabilities)
        prob_columns = target_labels['class_label'].unique()
        non_na_probs = target_labels['probability'].notna().sum()

        print(f"  Class labels found: {list(prob_columns)}")
        print(f"  Non-NaN probabilities: {non_na_probs}/{len(target_labels)}")

        # If most probabilities are NaN, this is likely a regression problem
        if non_na_probs < len(target_labels) * 0.1:  # Less than 10% valid probabilities
            print("  Detected regression problem - ROC curves not applicable")
            print(f"  Skipping '{target_value}' (use regression evaluation metrics instead)")
            continue

        # Debug: Check data quality
        total_rows = len(target_labels)
        missing_labels = target_labels['known_label'].isna().sum()
        missing_probs = target_labels['probability'].isna().sum()
        unique_samples = target_labels['sample_id'].nunique()
        unique_classes = target_labels['class_label'].nunique()

        print(f"  Data summary: {total_rows} total rows, {unique_samples} unique samples, {unique_classes} unique classes")
        print(f"  Missing data: {missing_labels} missing known_label, {missing_probs} missing probability")

        if missing_labels > 0:
            print(f"  Warning: Found {missing_labels} missing known_label values")
            missing_samples = target_labels[target_labels['known_label'].isna()]['sample_id'].unique()[:5]
            print(f"  Sample IDs with missing known_label: {list(missing_samples)}")

            # Remove rows with missing known_label
            target_labels = target_labels.dropna(subset=['known_label'])
            if target_labels.empty:
                print(f"  Error: No valid known_label data remaining for '{target_value}' - skipping")
                continue

        # 1. Pivot to wide format
        prob_df = target_labels.pivot(index='sample_id', columns='class_label', values='probability')

        print(f"  After pivot: {prob_df.shape[0]} samples x {prob_df.shape[1]} classes")
        print(f"  Class columns: {list(prob_df.columns)}")

        # Check for NaN values in probability data
        nan_counts = prob_df.isna().sum()
        if nan_counts.any():
            print(f"  NaN counts per class: {dict(nan_counts)}")
            print(f"  Samples with any NaN: {prob_df.isna().any(axis=1).sum()}/{len(prob_df)}")

            # Drop only rows where ALL probabilities are NaN
            all_nan_rows = prob_df.isna().all(axis=1)
            if all_nan_rows.any():
                print(f"  Dropping {all_nan_rows.sum()} samples with all NaN probabilities")
                prob_df = prob_df[~all_nan_rows]

            remaining_nans = prob_df.isna().sum().sum()
            if remaining_nans > 0:
                print(f"  Warning: {remaining_nans} individual NaN values remain - filling with 0")
                prob_df = prob_df.fillna(0)

            if prob_df.empty:
                print(f"  Error: No valid probability data remaining for '{target_value}' - skipping")
                continue

        # 2. Get true labels
        true_labels_df = target_labels.drop_duplicates('sample_id')[['sample_id', 'known_label']].set_index('sample_id')

        # 3. Align indices - only keep samples that exist in both datasets
        common_indices = prob_df.index.intersection(true_labels_df.index)
        if len(common_indices) == 0:
            print(f"  Error: No common sample_ids between probability and true label data for '{target_value}' - skipping")
            continue

        print(f"  Found {len(common_indices)} samples with both probability and true label data")

        # Filter both datasets to common indices
        prob_df_aligned = prob_df.loc[common_indices]
        y_true = true_labels_df.loc[common_indices]['known_label']

        # 4. Final check for NaN values
        if y_true.isna().any():
            print(f"  Error: True labels still contain NaN after alignment for '{target_value}' - skipping")
            continue

        if prob_df_aligned.isna().any().any():
            print(f"  Error: Probability data still contains NaN after alignment for '{target_value}' - skipping")
            continue

        # 5. Convert categorical labels to integer labels
        # Create a mapping from class names to integers
        class_names = list(prob_df_aligned.columns)
        class_to_int = {class_name: i for i, class_name in enumerate(class_names)}

        print(f"  Class mapping: {class_to_int}")

        # Convert true labels to integers
        y_true_np = y_true.map(class_to_int).to_numpy()
        y_probs_np = prob_df_aligned.to_numpy()

        print(f"  Data shape: y_true={y_true_np.shape}, y_probs={y_probs_np.shape}")
        print(f"  Unique true labels (integers): {set(y_true_np)}")
        print(f"  Class labels (columns): {class_names}")
        print(f"  Label distribution: {dict(zip(*np.unique(y_true_np, return_counts=True)))}")

        # Check for any unmapped labels (will be NaN)
        if pd.isna(y_true_np).any():
            print("  Error: Some true labels could not be mapped to class columns")
            unmapped_labels = set(y_true[y_true.map(class_to_int).isna()])
            print(f"  Unmapped labels: {unmapped_labels}")
            print(f"  Available classes: {class_names}")
            continue

        try:
            print(f"  Generating ROC curve for '{target_value}'...")
            fig = plot_roc_curves(y_true_np, y_probs_np)

            # Create output filename with target value
            safe_target_name = target_value.replace('/', '_').replace('\\', '_').replace(' ', '_')
            if len(target_values) > 1:
                output_filename = f"{output_name_base}_{safe_target_name}.{args.format}"
            else:
                output_filename = f"{output_name_base}.{args.format}"

            output_path = output_dir / output_filename
            print(f"  Saving ROC curve to: {output_path.absolute()}")
            fig.save(output_path, dpi=args.dpi, bbox_inches='tight')

        except Exception as e:
            print(f"  Error generating ROC curve for '{target_value}': {str(e)}")
            print(f"  Debug info - y_true type: {type(y_true_np)}, contains NaN: {pd.isna(y_true_np).any()}")
            print(f"  Debug info - y_probs type: {type(y_probs_np)}, contains NaN: {pd.isna(y_probs_np).any()}")
            continue

    print("ROC curves generated successfully!")


def generate_box_plots(labels, args, output_dir, output_name_base):
    """Generate box plots for model predictions"""

    print("Generating box plots...")

    # Parse target values from comma-separated string
    if args.target_value:
        target_values = [val.strip() for val in args.target_value.split(',')]
    else:
        # If no target values specified, use all unique variables
        target_values = labels['variable'].unique().tolist()

    print(f"Processing target values: {target_values}")

    for target_value in target_values:
        print(f"\nProcessing target value: '{target_value}'")

        # Filter labels for the current target value
        target_labels = labels[labels['variable'] == target_value]

        if target_labels.empty:
            print(f"  Warning: No data found for target value '{target_value}' - skipping")
            continue

        # Check if this is a classification problem (has probabilities)
        prob_columns = target_labels['class_label'].unique()
        non_na_probs = target_labels['probability'].notna().sum()

        print(f"  Class labels found: {list(prob_columns)}")
        print(f"  Non-NaN probabilities: {non_na_probs}/{len(target_labels)}")

        # If most probabilities are NaN, this is likely a regression problem
        if non_na_probs < len(target_labels) * 0.1:  # Less than 10% valid probabilities
            print("  Detected regression problem - precision-recall curves not applicable")
            print(f"  Skipping '{target_value}' (use regression evaluation metrics instead)")
            continue

        # Debug: Check data quality
        total_rows = len(target_labels)
        missing_labels = target_labels['known_label'].isna().sum()
        missing_probs = target_labels['probability'].isna().sum()
        unique_samples = target_labels['sample_id'].nunique()
        unique_classes = target_labels['class_label'].nunique()

        print(f"  Data summary: {total_rows} total rows, {unique_samples} unique samples, {unique_classes} unique classes")
        print(f"  Missing data: {missing_labels} missing known_label, {missing_probs} missing probability")

        if missing_labels > 0:
            print(f"  Warning: Found {missing_labels} missing known_label values")
            missing_samples = target_labels[target_labels['known_label'].isna()]['sample_id'].unique()[:5]
            print(f"  Sample IDs with missing known_label: {list(missing_samples)}")

            # Remove rows with missing known_label
            target_labels = target_labels.dropna(subset=['known_label'])
            if target_labels.empty:
                print(f"  Error: No valid known_label data remaining for '{target_value}' - skipping")
                continue

        # Remove rows with missing data
        clean_data = target_labels.dropna(subset=['known_label', 'probability'])

        if clean_data.empty:
            print("    No valid data after cleaning - skipping")
            continue

        # Get unique classes
        classes = clean_data['class_label'].unique()

        for class_label in classes:
            print(f"    Generating box plot for class: {class_label}")

            # Filter for current class
            class_data = clean_data[clean_data['class_label'] == class_label]

            try:
                # Create the box plot
                fig = plot_boxplot(
                    categorical_x=class_data['known_label'],
                    numerical_y=class_data['probability'],
                    title_x='True Label',
                    title_y=f'Predicted Probability ({class_label})',
                )

                # Save the plot
                safe_class_name = str(class_label).replace('/', '_').replace('\\', '_').replace(' ', '_').replace(':', '_')
                safe_target_name = target_value.replace('/', '_').replace('\\', '_').replace(' ', '_')
                output_filename = f"{output_name_base}_{safe_target_name}_{safe_class_name}.{args.format}"
                output_path = output_dir / output_filename

                print(f"      Saving box plot to: {output_path.absolute()}")
                fig.savefig(output_path, dpi=args.dpi, bbox_inches='tight')
                plt.close(fig)

            except Exception as e:
                print(f"      Error generating box plot for class '{class_label}': {str(e)}")
                continue


def main():
    """Main function to parse arguments and generate plots"""
    parser = argparse.ArgumentParser(description="Generate plots using flexynesis")

    # Required arguments
    parser.add_argument("--labels", type=str, required=False,
                        help="Path to labels file generated by flexynesis")

    # Plot type
    parser.add_argument("--plot_type", type=str, required=True,
                        choices=['dimred', 'kaplan_meier', 'cox', 'scatter', 'concordance_heatmap', 'pr_curve', 'roc_curve', 'box_plot'],
                        help="Type of plot to generate: 'dimred' for dimensionality reduction, 'kaplan_meier' for survival analysis, 'cox' for Cox proportional hazards analysis, 'scatter' for scatter plots, 'concordance_heatmap' for label concordance heatmaps, 'pr_curve' for precision-recall curves, 'roc_curve' for ROC curves, or 'box_plot' for box plots.")

    # Arguments for dimensionality reduction
    parser.add_argument("--embeddings", type=str,
                        help="Path to input data embeddings file (CSV or tabular format). Required for dimred plots.")
    parser.add_argument("--method", type=str, default='pca', choices=['pca', 'umap'],
                        help="Transformation method ('pca' or 'umap'). Default is 'pca'. Used for dimred plots.")

    # Arguments for Kaplan-Meier
    parser.add_argument("--survival_data", type=str,
                        help="Path to survival data file with columns: duration and event. Required for kaplan_meier plots.")
    parser.add_argument("--surv_time_var", type=str, required=False,
                        help="Column name for survival time")
    parser.add_argument("--surv_event_var", type=str, required=False,
                        help="Column name for survival event")

    # Arguments for Cox analysis
    parser.add_argument("--important_features", type=str,
                        help="Path to calculated feature importance file. Required for cox plots.")
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
    parser.add_argument("--crossval", action='store_true',
                        help="If True, performs K-fold cross-validation and returns average C-index. Default is False")
    parser.add_argument("--n_splits", type=int, default=5,
                        help="Number of folds for cross-validation. Default is 5")
    parser.add_argument("--random_state", type=int, default=42,
                        help="Random seed for reproducibility. Default is 42")

    # Arguments for dimred, scatter plot, heatmap, PR curves, ROC curves, and box plots
    parser.add_argument("--target_value", type=str, default=None,
                        help="Target value for scatter plot.")

    # Arguments for dimred
    parser.add_argument("--color", type=str, default=None,
                        help="User-defined color for the plot.")

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
        if args.plot_type not in ['dimred', 'kaplan_meier', 'cox', 'scatter', 'concordance_heatmap', 'pr_curve', 'roc_curve', 'box_plot']:
            raise ValueError(f"Invalid plot type: {args.plot_type}. Must be one of: 'dimred', 'kaplan_meier', 'cox', 'scatter', 'concordance_heatmap', 'pr_curve', 'roc_curve', 'box_plot'")

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

        if args.plot_type in ['cox']:
            if not args.important_features:
                raise ValueError("--important_features is required when plot_type is 'cox'")
            if not os.path.isfile(args.important_features):
                raise FileNotFoundError(f"Important features file not found: {args.important_features}")
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
            if not args.crossval:
                args.crossval = False
            if not isinstance(args.n_splits, int) or args.n_splits <= 0:
                raise ValueError("--n_splits must be a positive integer")
            if not isinstance(args.random_state, int):
                raise ValueError("--random_state must be an integer")

        if args.plot_type in ['scatter']:
            if not args.labels:
                raise ValueError("--labels is required for scatter plots")
            if not args.target_value:
                print("--target_value is not specified, using all unique variables from labels")
            if not os.path.isfile(args.labels):
                raise FileNotFoundError(f"Labels file not found: {args.labels}")

        if args.plot_type in ['concordance_heatmap']:
            if not args.labels:
                raise ValueError("--labels is required for concordance heatmap")
            if not args.target_value:
                print("--target_value is not specified, using all unique variables from labels")
            if not os.path.isfile(args.labels):
                raise FileNotFoundError(f"Labels file not found: {args.labels}")

        if args.plot_type in ['pr_curve']:
            if not args.labels:
                raise ValueError("--labels is required for precision-recall curves")
            if not args.target_value:
                print("--target_value is not specified, using all unique variables from labels")
            if not os.path.isfile(args.labels):
                raise FileNotFoundError(f"Labels file not found: {args.labels}")

        if args.plot_type in ['roc_curve']:
            if not args.labels:
                raise ValueError("--labels is required for ROC curves")
            if not args.target_value:
                print("--target_value is not specified, using all unique variables from labels")
            if not os.path.isfile(args.labels):
                raise FileNotFoundError(f"Labels file not found: {args.labels}")

        if args.plot_type in ['box_plot']:
            if not args.labels:
                raise ValueError("--labels is required for box plots")
            if not args.target_value:
                print("--target_value is not specified, using all unique variables from labels")
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
            elif args.plot_type == 'pr_curve':
                labels_name = Path(args.labels).stem
                output_name_base = f"{labels_name}_pr_curves"
            elif args.plot_type == 'roc_curve':
                labels_name = Path(args.labels).stem
                output_name_base = f"{labels_name}_roc_curves"
            elif args.plot_type == 'box_plot':
                labels_name = Path(args.labels).stem
                output_name_base = f"{labels_name}_box_plot"

        # Generate plots based on type
        if args.plot_type in ['dimred']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)
            # Load embeddings data
            print(f"Loading embeddings from: {args.embeddings}")
            embeddings, sample_names = load_embeddings(args.embeddings)
            print(f"embeddings shape: {embeddings.shape}")

            # Match samples to embeddings
            matched_labels = match_samples_to_embeddings(sample_names, label_data)
            print(f"Successfully matched {len(matched_labels)} samples for dimensionality reduction")
            print(f"Matched labels shape: {matched_labels.shape}")
            print(f"Columns in matched labels: {matched_labels.columns.tolist()}")
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
            # Load important_features and datasets
            print(f"Loading important features from: {args.important_features}")
            important_features = load_labels(args.important_features)
            print(f"Loading training dataset from: {args.clinical_train}")
            clinical_train = load_omics(args.clinical_train)
            print(f"Loading test dataset from: {args.clinical_test}")
            clinical_test = load_omics(args.clinical_test)
            print(f"Loading training omics dataset from: {args.omics_train}")
            omics_train = load_omics(args.omics_train)
            print(f"Loading test omics dataset from: {args.omics_test}")
            omics_test = load_omics(args.omics_test)

            generate_cox_plots(important_features, clinical_train, clinical_test, omics_test, omics_train, args, output_dir, output_name_base)

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

        elif args.plot_type in ['pr_curve']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)

            generate_pr_curves(label_data, args, output_dir, output_name_base)

        elif args.plot_type in ['roc_curve']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)

            generate_roc_curves(label_data, args, output_dir, output_name_base)

        elif args.plot_type in ['box_plot']:
            # Load labels
            print(f"Loading labels from: {args.labels}")
            label_data = load_labels(args.labels)

            generate_box_plots(label_data, args, output_dir, output_name_base)

        print("All plots generated successfully!")

    except (FileNotFoundError, ValueError, pd.errors.ParserError) as e:
        print(f"Error: {e}")
        return 1

    return 0


if __name__ == "__main__":
    exit(main())
