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
    def __init__(self, dataset, task, n_splits=3, quality_threshold=0.2, knn_k=10):
        self.dataset = dataset
        self.task = task
        self.n_splits = n_splits
        self.quality_threshold = quality_threshold
        self.knn_k = knn_k
        self.issues = None
        self.knn_graph = None
        self.features = self.dataset.drop('target', axis=1).columns.tolist()
        self.issue_summary = None

    def report_issues(self):
        X = self.dataset.drop('target', axis=1)
        y = self.dataset['target']

        # Ensure compatibility with Galaxy
        X = X.to_numpy() if hasattr(X, 'to_numpy') else np.asarray(X)
        y = y.to_numpy() if hasattr(y, 'to_numpy') else np.asarray(y)

        if self.task == 'classification':
            model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
            cv = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=42)
            pred_probs = cross_val_predict(model, X, y, cv=cv, method='predict_proba')

            lab = Datalab(self.dataset, label_name='target')
            lab.find_issues(pred_probs=pred_probs)
            self.issues = lab.get_issues()
            self.issue_summary = lab.get_issue_summary()
            print(self.issue_summary)

        elif self.task == 'regression':
            model = LinearRegression()
            cv = KFold(n_splits=self.n_splits, shuffle=True, random_state=42)
            pred_y = cross_val_predict(model, X, y, cv=cv, method='predict')
            scores = get_label_quality_scores(y, pred_y, method='residual')
            self.issues = pd.DataFrame({
                'label_quality': scores,
                'is_label_issue': scores < self.quality_threshold
            })

        return self.dataset.copy(), self.issues.copy()

    def clean_selected_issues(self, method='remove', label_issues=True, outliers=True, near_duplicates=True, non_iid=True):
        if self.issues is None:
            raise RuntimeError("Must run report_issues() before cleaning.")

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
            fixed.loc[to_fix, 'target'] = most_likely[to_fix]
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
    parser.add_argument("--method", default="remove", choices=["remove", "replace"], help="Cleaning method")
    parser.add_argument("--summary", action="store_true", help="Print and save issue summary only, no cleaning")
    parser.add_argument("--no-label-issues", action="store_true", help="Exclude label issues from cleaning")
    parser.add_argument("--no-outliers", action="store_true", help="Exclude outlier issues from cleaning")
    parser.add_argument("--no-near-duplicates", action="store_true", help="Exclude near-duplicate issues from cleaning")
    parser.add_argument("--no-non-iid", action="store_true", help="Exclude non-i.i.d. issues from cleaning")

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
    handler = IssueHandler(dataset=df, task=args.task)
    _, issues = handler.report_issues()

    # Save summary
    if handler.issue_summary is not None:
        with open("summary.txt", "w") as f:
            f.write(str(handler.issue_summary))
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
