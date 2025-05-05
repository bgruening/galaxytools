import argparse

import numpy as np
import pandas as pd
from xgboost import XGBClassifier

from cleanlab import Datalab
from cleanlab.regression.rank import get_label_quality_scores
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import KFold, StratifiedKFold, cross_val_predict
from sklearn.neighbors import NearestNeighbors


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
        self.features = None
        self.knn_graph = None
        self.features = self.dataset.drop('target', axis=1).columns.tolist()
        self.pred_probs = None
        self.issue_summary = None

    def report_issues(self):
        X = self.dataset.drop('target', axis=1)
        y = self.dataset['target']

        # Ensure compatibility with Galaxy!
        X = X.to_numpy() if hasattr(X, 'to_numpy') else np.asarray(X)
        y = y.to_numpy() if hasattr(y, 'to_numpy') else np.asarray(y)

        # Compute knn_graph
        nn = NearestNeighbors(n_neighbors=self.knn_k + 1)
        nn.fit(X)
        self.knn_graph = nn.kneighbors(return_distance=False)[:, 1:]

        if self.task == 'classification':
            model = XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
            cv = StratifiedKFold(n_splits=self.n_splits, shuffle=True, random_state=42)
            self.pred_probs = cross_val_predict(model, X, y, cv=cv, method='predict_proba')

            lab = Datalab(self.dataset, label_name='target')
            lab.find_issues(pred_probs=self.pred_probs)  # simplified for compatibility. in next version, knn_graph and features are going to be used to acquire more issue types.
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

    def clean_selected_issues(
        self,
        method='remove',
        label_issues=True,
        outliers=True,
        near_duplicates=True,
        non_iid=True
    ):
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
    parser.add_argument("--csv", required=True, help="Path to dataset CSV (must include a 'target' column)")
    parser.add_argument("--task", required=True, choices=["classification", "regression"], help="Type of ML task")
    parser.add_argument("--method", default="remove", choices=["remove", "replace"], help="Cleaning method")
    parser.add_argument("--summary", action="store_true", help="Print and save issue summary only, no cleaning")

    parser.add_argument("--no-label-issues", action="store_true", help="Exclude label issues from cleaning")
    parser.add_argument("--no-outliers", action="store_true", help="Exclude outlier issues from cleaning")
    parser.add_argument("--no-near-duplicates", action="store_true", help="Exclude near-duplicate issues from cleaning")
    parser.add_argument("--no-non-iid", action="store_true", help="Exclude non-i.i.d. issues from cleaning")

    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    if 'target' not in df.columns:
        raise ValueError("Dataset must contain a 'target' column.")

    handler = IssueHandler(dataset=df, task=args.task)
    _, issues = handler.report_issues()

    if handler.issue_summary is not None:
        with open("summary.txt", "w") as f:
            f.write(str(handler.issue_summary))
        print("Issue summary saved to: summary.txt")

    if args.summary:
        return

    cleaned_df = handler.clean_selected_issues(
        method=args.method,
        label_issues=not args.no_label_issues,
        outliers=not args.no_outliers,
        near_duplicates=not args.no_near_duplicates,
        non_iid=not args.no_non_iid
    )

    cleaned_df.to_csv("cleaned_data.csv", index=False)
    print("Cleaned dataset saved to: cleaned_data.csv")


# -------------------
# Entry point
# -------------------
if __name__ == "__main__":
    main()
