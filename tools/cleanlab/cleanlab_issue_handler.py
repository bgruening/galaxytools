import argparse

import numpy as np
import pandas as pd
from cleanlab.datalab.datalab import Datalab
from cleanlab.regression.rank import get_label_quality_scores
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_predict, KFold, StratifiedKFold
from xgboost import XGBClassifier

# -------------------
# Issue Handler
# -------------------


class IssueHandler:
    def __init__(self, dataset, task, target_column, n_splits=3, quality_threshold=0.2, knn_k=10):
        self.dataset = dataset
        self.task = task
        self.target_column = target_column
        self.n_splits = n_splits
        self.quality_threshold = quality_threshold
        self.knn_k = knn_k
        self.issues = None
        self.knn_graph = None
        self.features = self.dataset.drop(target_column, axis=1).columns.tolist()
        self.issue_summary = None
        self.pred_probs = None

    def report_issues(self):
        X = self.dataset.drop(self.target_column, axis=1)
        y = self.dataset[self.target_column]

        # Ensure compatibility with Galaxy
        X = X.to_numpy() if hasattr(X, 'to_numpy') else np.asarray(X)
        y = y.to_numpy() if hasattr(y, 'to_numpy') else np.asarray(y)

        if self.task == 'classification':
            model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
            cv = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=42)
            self.pred_probs = cross_val_predict(model, X, y, cv=cv, method='predict_proba')

            lab = Datalab(self.dataset, label_name=self.target_column)
            lab.find_issues(pred_probs=self.pred_probs)
            self.issues = lab.get_issues()
            self.issue_summary = lab.get_issue_summary()
            print(self.issue_summary)

        elif self.task == 'regression':
            model = LinearRegression()
            cv = KFold(n_splits=self.n_splits, shuffle=True, random_state=42)
            pred_y = cross_val_predict(model, X, y, cv=cv, method='predict')
            scores = get_label_quality_scores(y, pred_y, method='residual')
            is_low_quality = scores < self.quality_threshold
            self.issues = pd.DataFrame({
                'label_quality': scores,
                'is_low_quality': is_low_quality
            })
            self.issue_summary = {
                'quality_threshold': self.quality_threshold,
                'num_low_quality': int(is_low_quality.sum()),
                'mean_label_quality': float(np.mean(scores)),
                'median_label_quality': float(np.median(scores)),
                'min_label_quality': float(np.min(scores)),
                'max_label_quality': float(np.max(scores)),
            }
            print("Regression Issue Summary:")
            for k, v in self.issue_summary.items():
                print(f"{k.replace('_', ' ').capitalize()}: {v:.4f}" if isinstance(v, float) else f"{k.replace('_', ' ').capitalize()}: {v}")

        return self.dataset.copy(), self.issues.copy(), self.issue_summary

    def clean_selected_issues(self, method='remove', label_issues=True, outliers=True, near_duplicates=True, non_iid=True):
        if self.issues is None:
            raise RuntimeError("Must run report_issues() before cleaning.")

        if self.task == 'regression':
            clean_mask = self.issues['is_low_quality'].fillna(False)
        else:
            clean_mask = pd.Series([False] * len(self.dataset))
            for issue_type, use_flag in [
                ('is_label_issue', label_issues),
                ('is_outlier_issue', outliers),
                ('is_near_duplicate_issue', near_duplicates),
                ('is_non_iid_issue', non_iid)
            ]:
                if use_flag and issue_type in self.issues.columns:
                    clean_mask |= self.issues[issue_type].fillna(False)

        if method == 'remove':
            return self.dataset[~clean_mask].copy()

        elif method == 'replace' and self.task == 'classification':
            most_likely = np.argmax(self.pred_probs, axis=1)
            fixed = self.dataset.copy()
            to_fix = self.issues['is_label_issue'] & label_issues
            fixed.loc[to_fix, self.target_column] = most_likely[to_fix]
            return fixed

        elif method == 'replace' and self.task == 'regression':
            raise NotImplementedError("Replace method not implemented for regression label correction.")

        else:
            raise ValueError("Invalid method or unsupported combination.")

# -------------------
# Main CLI Entry
# -------------------

def main():
    parser = argparse.ArgumentParser(description="Cleanlab Issue Handler CLI")
    parser.add_argument("--input_file", nargs=2, required=True, metavar=('FILE', 'EXT'), help="Input file path and its extension")
    parser.add_argument("--task", required=True, choices=["classification", "regression"], help="Type of ML task")
    parser.add_argument("--target_column", default="target", help="Name of the target column")
    parser.add_argument("--method", default="remove", choices=["remove", "replace"], help="Cleaning method")
    parser.add_argument("--summary", action="store_true", help="Print and save issue summary only, no cleaning")
    parser.add_argument("--no-label-issues", action="store_true", help="Exclude label issues from cleaning")
    parser.add_argument("--no-outliers", action="store_true", help="Exclude outlier issues from cleaning")
    parser.add_argument("--no-near-duplicates", action="store_true", help="Exclude near-duplicate issues from cleaning")
    parser.add_argument("--no-non-iid", action="store_true", help="Exclude non-i.i.d. issues from cleaning")
    parser.add_argument('--quality-threshold', type=float, default=0.2, help='Threshold for low-quality labels (regression only)')

    args = parser.parse_args()

    # Load dataset based on file extension
    file_path, file_ext = args.input_file
    file_ext = file_ext.lower()

    print(f"Loading dataset from: {file_path} with extension: {file_ext}")

    if file_ext == "csv":
        df = pd.read_csv(file_path)
    elif file_ext in ["tsv", "tabular"]:
        df = pd.read_csv(file_path, sep="\t")
    else:
        raise ValueError(f"Unsupported file format: {file_ext}")

    # Run IssueHandler
    handler = IssueHandler(dataset=df, 
                           task=args.task, 
                           target_column=args.target_column, 
                           quality_threshold=args.quality_threshold)
    _, issues, summary = handler.report_issues()

    # Save summary
    if summary is not None:
        with open("summary.txt", "w") as f:
            if args.task == "regression":
                f.write("Regression Issue Summary:\n")
                for k, v in summary.items():
                    text = f"{k.replace('_', ' ').capitalize()}: {v:.4f}" if isinstance(v, float) else f"{k.replace('_', ' ').capitalize()}: {v}"
                    f.write(text + "\n")
            else:
                f.write(str(summary))
        print("Issue summary saved to: summary.txt")

    if args.summary:
        return

    # Clean selected issues
    cleaned_df = handler.clean_selected_issues(
        method=args.method,
        label_issues=not args.no_label_issues,
        outliers=not args.no_outliers,
        near_duplicates=not args.no_near_duplicates,
        non_iid=not args.no_non_iid
    )

    print(f"Cleaned dataset shape: {cleaned_df.shape}")
    print(f"Original dataset shape: {df.shape}")

    output_filename = "cleaned_data"
    if file_ext == "csv":
        cleaned_df.to_csv(output_filename, index=False)
    elif file_ext in ["tsv", "tabular"]:
        cleaned_df.to_csv(output_filename, sep="\t", index=False)
    else:
        raise ValueError(f"Unsupported output format: {file_ext}")

    print(f"Cleaned dataset saved to: {output_filename}")

# -------------------
# Entry point
# -------------------


if __name__ == "__main__":
    main()
