"""
Tabular data prediction using TabPFN
"""
import argparse
import time

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import (
    average_precision_score,
    precision_recall_curve,
    root_mean_squared_error,
    r2_score,
)
from tabpfn import TabPFNClassifier, TabPFNRegressor


def separate_features_labels(data):
    df = pd.read_csv(data, sep="\t")
    labels = df.iloc[:, -1]
    features = df.iloc[:, :-1]
    return features, labels


def classification_plot(xval, yval, leg_label, title, xlabel, ylabel):
    plt.figure(figsize=(8, 6))
    plt.plot(xval, yval, label=leg_label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
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
    tr_features, tr_labels = separate_features_labels(args["train_data"])
    # prepare test data
    if args["testhaslabels"] == "haslabels":
        te_features, te_labels = separate_features_labels(args["test_data"])
    else:
        te_features = pd.read_csv(args["test_data"], sep="\t")
        te_labels = []
    s_time = time.time()
    if args["selected_task"] == "Classification":
        classifier = TabPFNClassifier(device="cpu")
        classifier.fit(tr_features, tr_labels)
        y_eval = classifier.predict(te_features)
        pred_probas_test = classifier.predict_proba(te_features)
        if len(te_labels) > 0:
            precision, recall, thresholds = precision_recall_curve(
                te_labels, pred_probas_test[:, 1]
            )
            average_precision = average_precision_score(
                te_labels, pred_probas_test[:, 1]
            )
            classification_plot(
                recall,
                precision,
                f"Precision-Recall Curve (AP={average_precision:.2f})",
                "Precision-Recall Curve",
                "Recall",
                "Precision",
            )
    else:
        regressor = TabPFNRegressor(device="cpu")
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
    e_time = time.time()
    print(
        "Time taken by TabPFN for training and prediction: {} seconds".format(
            e_time - s_time
        )
    )
    te_features["predicted_labels"] = y_eval
    te_features.to_csv("output_predicted_data", sep="\t", index=None)


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
