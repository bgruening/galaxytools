"""
Tabular data prediction using TabPFN
"""
import argparse
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from catboost import CatBoostClassifier, CatBoostRegressor
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    r2_score,
    root_mean_squared_error,
)
from sklearn.preprocessing import label_binarize
from tabpfn import TabPFNClassifier, TabPFNRegressor


def separate_features_labels(data):
    df = pd.read_csv(data, sep="\t")
    labels = df.iloc[:, -1]
    features = df.iloc[:, :-1]
    return features, labels


def classification_plot(y_true, y_scores, m_name):
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
        plt.title(f"{m_name}: Precision-Recall Curve (binary classification)")
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
        plt.title(
            "{}: Precision-Recall Curve (Multiclass Classification)".format(m_name)
        )
    plt.xlabel("Recall")
    plt.ylabel("Precision")
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.savefig(f"output_plot_{m_name}.png")


def regression_plot(xval, yval, title, xlabel, ylabel, m_name):
    plt.figure(figsize=(8, 6))
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend(loc="lower left")
    plt.grid(True)
    plt.scatter(xval, yval, alpha=0.8)
    xticks = np.arange(len(xval))
    plt.plot(xticks, xticks, color="red", linestyle="--", label="y = x")
    plt.savefig("output_plot_{}.png".format(m_name))


def train_evaluate(args):
    """
    Train TabPFN and predict
    """
    # prepare train data
    tr_features, tr_labels = separate_features_labels(args["train_data"])
    # prepare test data
    if args["testhaslabels"] == "true":
        te_features, te_labels = separate_features_labels(args["test_data"])
    else:
        te_features = pd.read_csv(args["test_data"], sep="\t")
        te_labels = []
    s_time = time.time()
    if args["selected_task"] == "Classification":
        models = [
            ("TabPFN", TabPFNClassifier(random_state=42)),
            ("CatBoost", CatBoostClassifier(random_state=42, verbose=0)),
        ]
        for m_name, model in models:
            model.fit(tr_features, tr_labels)
            y_eval = model.predict(te_features)
            pred_probas_test = model.predict_proba(te_features)
            if len(te_labels) > 0:
                classification_plot(te_labels, pred_probas_test, m_name)
            te_features["predicted_labels"] = y_eval
            te_features.to_csv(
                "output_predicted_data_{}".format(m_name), sep="\t", index=None
            )
    else:
        models = [
            ("TabPFN", TabPFNRegressor(random_state=42)),
            ("CatBoost", CatBoostRegressor(random_state=42, verbose=0)),
        ]
        for m_name, model in models:
            model.fit(tr_features, tr_labels)
            y_eval = model.predict(te_features)
            if len(te_labels) > 0:
                score = root_mean_squared_error(te_labels, y_eval)
                r2_metric_score = r2_score(te_labels, y_eval)
                regression_plot(
                    te_labels,
                    y_eval,
                    f"Scatter plot for {m_name}: True vs predicted values. RMSE={score:.2f}, R2={r2_metric_score:.2f}",
                    "True values",
                    "Predicted values",
                    m_name,
                )
            te_features["predicted_labels"] = y_eval
            te_features.to_csv(
                "output_predicted_data_{}".format(m_name), sep="\t", index=None
            )
    e_time = time.time()
    print(
        "Time taken by TabPFN for training and prediction: {} seconds".format(
            e_time - s_time
        )
    )


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-trdata", "--train_data", required=True, help="Train data")
    arg_parser.add_argument("-tedata", "--test_data", required=True, help="Test data")
    arg_parser.add_argument(
        "-testhaslabels",
        "--testhaslabels",
        required=True,
        help="if test data contain labels",
    )
    arg_parser.add_argument(
        "-selectedtask",
        "--selected_task",
        required=True,
        help="Type of machine learning task",
    )
    # get argument values
    args = vars(arg_parser.parse_args())
    train_evaluate(args)
