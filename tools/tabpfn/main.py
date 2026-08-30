"""
Tabular data prediction using TabPFN
"""

import argparse
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    accuracy_score,
    average_precision_score,
    balanced_accuracy_score,
    f1_score,
    precision_recall_curve,
    r2_score,
    root_mean_squared_error,
)
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.preprocessing import label_binarize
from tabpfn import TabPFNClassifier, TabPFNRegressor

MAX_IGNORE_PRETRAINING_LIMITS_SAMPLES = 1000
SEED = 42


def separate_features_labels(data, header):
    df = pd.read_csv(data, sep="\t", header=0 if header == "true" else None)
    labels = df.iloc[:, -1]
    features = df.iloc[:, :-1]
    return features, labels


def make_estimator(selected_task, model_path, n_samples):
    estimator_class = (
        TabPFNClassifier if selected_task == "Classification" else TabPFNRegressor
    )
    kwargs = {
        "random_state": SEED,
        "model_path": model_path,
    }
    if n_samples > MAX_IGNORE_PRETRAINING_LIMITS_SAMPLES:
        kwargs["ignore_pretraining_limits"] = True
    return estimator_class(**kwargs)


def classification_plot(y_true, y_scores):
    plt.figure(figsize=(8, 6))
    is_binary = len(np.unique(y_true)) == 2
    if is_binary:
        # Compute precision-recall curve
        precision, recall, _ = precision_recall_curve(y_true, y_scores[:, 1])
        average_precision = average_precision_score(y_true, y_scores[:, 1])
        plt.plot(
            recall,
            precision,
            label=f"Precision-Recall Curve (AP={average_precision:.2f})",
        )
        plt.title("Precision-Recall Curve (binary classification)")
    else:
        y_true_bin = label_binarize(y_true, classes=np.unique(y_true))
        n_classes = y_true_bin.shape[1]
        class_labels = [f"Class {i}" for i in range(n_classes)]
        # Plot PR curve for each class
        for i in range(n_classes):
            precision, recall, _ = precision_recall_curve(
                y_true_bin[:, i], y_scores[:, i]
            )
            ap_score = average_precision_score(y_true_bin[:, i], y_scores[:, i])
            plt.plot(
                recall, precision, label=f"{class_labels[i]} (AP = {ap_score:.2f})"
            )
        # Compute micro-average PR curve
        precision, recall, _ = precision_recall_curve(
            y_true_bin.ravel(), y_scores.ravel()
        )
        plt.plot(
            recall, precision, linestyle="--", color="black", label="Micro-average"
        )
        plt.title("Precision-Recall Curve (Multiclass Classification)")
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.savefig("output_plot.png")


def regression_plot(xval, yval, title, xlabel, ylabel):
    plt.figure(figsize=(8, 6))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.scatter(xval, yval, alpha=0.8)
    xticks = np.arange(len(xval))
    plt.plot(xticks, xticks, color="red", linestyle="--", label="y = x")
    plt.savefig("output_plot.png")


def train_evaluate(args):
    """
    Train TabPFN and predict
    """
    # prepare train data
    tr_features, tr_labels = separate_features_labels(
        args["train_data"], args["train_header"]
    )
    # prepare test data
    if args["testhaslabels"] == "true":
        te_features, te_labels = separate_features_labels(
            args["test_data"], args["test_header"]
        )
    else:
        te_features = pd.read_csv(
            args["test_data"],
            sep="\t",
            header=0 if args["test_header"] == "true" else None,
        )
        te_labels = []
    s_time = time.time()
    if args["selected_task"] == "Classification":
        classifier = make_estimator(
            args["selected_task"], args["model_path"], tr_features.shape[0]
        )
        classifier.fit(tr_features, tr_labels)
        y_eval = classifier.predict(te_features)
        pred_probas_test = classifier.predict_proba(te_features)
        if len(te_labels) > 0:
            classification_plot(te_labels, pred_probas_test)
        te_features["predicted_labels"] = y_eval
        te_features.to_csv("output_predicted_data", sep="\t", index=None)
    else:
        regressor = make_estimator(
            args["selected_task"], args["model_path"], tr_features.shape[0]
        )
        regressor.fit(tr_features, tr_labels)
        y_eval = regressor.predict(te_features)
        if len(te_labels) > 0:
            score = root_mean_squared_error(te_labels, y_eval)
            r2_metric_score = r2_score(te_labels, y_eval)
            regression_plot(
                te_labels,
                y_eval,
                f"Scatter plot: True vs predicted values. RMSE={score:.2f}, R2={r2_metric_score:.2f}",
                "True values",
                "Predicted values",
            )
    te_features["predicted_labels"] = y_eval
    te_features.to_csv("output_predicted_data", sep="\t", index=None)
    e_time = time.time()
    print(
        f"Time taken by TabPFN for training and prediction: {e_time - s_time} seconds"
    )


def append_summary_rows(metrics_df, metric_columns):
    summary_rows = []
    for summary_name, summary_function in (("mean", "mean"), ("std", "std")):
        row = {"fold": summary_name}
        for metric_column in metric_columns:
            row[metric_column] = getattr(metrics_df[metric_column], summary_function)()
        summary_rows.append(row)
    return pd.concat([metrics_df, pd.DataFrame(summary_rows)], ignore_index=True)


def cross_validate(args):
    """
    Evaluate TabPFN via cross-validation and export predictions and metrics.
    """
    n_splits = int(args["n_splits"])
    cv_strategy = args["cv_strategy"]
    features, labels = separate_features_labels(
        args["train_data"], args["train_header"]
    )
    predictions = pd.Series(index=features.index, dtype=object)
    fold_numbers = pd.Series(index=features.index, dtype="Int64")
    metrics = []
    s_time = time.time()

    if args["selected_task"] == "Classification":
        if cv_strategy == "stratified":
            class_counts = labels.value_counts()
            too_small_classes = class_counts[class_counts < n_splits]
            if not too_small_classes.empty:
                too_small_class_names = ", ".join(
                    str(name) for name in too_small_classes.index
                )
                raise ValueError(
                    "Cannot run stratified cross validation: each class must contain at "
                    f"least {n_splits} samples. Classes below that limit: "
                    f"{too_small_class_names}"
                )
            splitter = StratifiedKFold(
                n_splits=n_splits, shuffle=True, random_state=SEED
            )
            split_iterator = splitter.split(features, labels)
        else:
            splitter = KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
            split_iterator = splitter.split(features)
        for fold_index, (train_index, test_index) in enumerate(split_iterator, start=1):
            estimator = make_estimator(
                args["selected_task"], args["model_path"], len(train_index)
            )
            estimator.fit(features.iloc[train_index], labels.iloc[train_index])
            fold_predictions = estimator.predict(features.iloc[test_index])
            predictions.iloc[test_index] = fold_predictions
            fold_numbers.iloc[test_index] = fold_index
            metrics.append(
                {
                    "fold": fold_index,
                    "accuracy": accuracy_score(
                        labels.iloc[test_index], fold_predictions
                    ),
                    "balanced_accuracy": balanced_accuracy_score(
                        labels.iloc[test_index], fold_predictions
                    ),
                    "f1_weighted": f1_score(
                        labels.iloc[test_index],
                        fold_predictions,
                        average="weighted",
                        zero_division=0,
                    ),
                }
            )
        metric_columns = ["accuracy", "balanced_accuracy", "f1_weighted"]
    else:
        if cv_strategy == "stratified":
            raise ValueError(
                "Stratified cross validation is only available for classification. "
                "Use kfold cross validation for regression."
            )
        splitter = KFold(n_splits=n_splits, shuffle=True, random_state=SEED)
        for fold_index, (train_index, test_index) in enumerate(
            splitter.split(features), start=1
        ):
            estimator = make_estimator(
                args["selected_task"], args["model_path"], len(train_index)
            )
            estimator.fit(features.iloc[train_index], labels.iloc[train_index])
            fold_predictions = estimator.predict(features.iloc[test_index])
            predictions.iloc[test_index] = fold_predictions
            fold_numbers.iloc[test_index] = fold_index
            metrics.append(
                {
                    "fold": fold_index,
                    "rmse": root_mean_squared_error(
                        labels.iloc[test_index], fold_predictions
                    ),
                    "r2": r2_score(labels.iloc[test_index], fold_predictions),
                }
            )
        metric_columns = ["rmse", "r2"]

    output_data = features.copy()
    output_data["true_labels"] = labels
    output_data["fold"] = fold_numbers
    output_data["predicted_labels"] = predictions
    output_data.to_csv("output_predicted_data", sep="\t", index=None)

    metrics_df = pd.DataFrame(metrics)
    metrics_df = append_summary_rows(metrics_df, metric_columns)
    metrics_df.to_csv("cv_metrics.tsv", sep="\t", index=None)

    e_time = time.time()
    print(f"Time taken by TabPFN for cross validation: {e_time - s_time} seconds")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument(
        "-evalmethod",
        "--evaluation_method",
        required=True,
        choices=["train_test", "cross_validation"],
        help="Evaluation method",
    )
    arg_parser.add_argument("-trdata", "--train_data", required=True, help="Train data")
    arg_parser.add_argument("-tedata", "--test_data", help="Test data")
    arg_parser.add_argument(
        "-train_header",
        "--train_header",
        required=True,
        help="if train data contain header",
    )
    arg_parser.add_argument(
        "-test_header", "--test_header", help="if test data contain header"
    )
    arg_parser.add_argument(
        "-testhaslabels",
        "--testhaslabels",
        help="if test data contain labels",
    )
    arg_parser.add_argument(
        "-nsplits",
        "--n_splits",
        type=int,
        default=5,
        help="Number of cross-validation folds",
    )
    arg_parser.add_argument(
        "-cvstrategy",
        "--cv_strategy",
        choices=["kfold", "stratified"],
        default="stratified",
        help="Cross-validation splitting strategy",
    )
    arg_parser.add_argument(
        "-selectedtask",
        "--selected_task",
        required=True,
        help="Type of machine learning task",
    )
    arg_parser.add_argument(
        "-modelpath",
        "--model_path",
        required=True,
        help="Pretrained model to use",
    )
    # get argument values
    args = vars(arg_parser.parse_args())
    if args["evaluation_method"] == "cross_validation":
        cross_validate(args)
    else:
        train_evaluate(args)
