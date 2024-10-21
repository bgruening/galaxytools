"""
Tabular data prediction using TabPFN
"""

import argparse
import time

from sklearn.metrics import accuracy_score
from tabpfn import TabPFNClassifier
import pandas as pd
import torch


def separate_features_labels(data):
    df = pd.read_csv(data, sep=",")
    labels = df.iloc[:, -1]
    features = df.iloc[:, :-1]
    print(df)
    print(features)
    print(labels)
    return features, labels


def train_evaluate(args):
    """
    Train TabPFN
    """
    print(args)

    tr_features, tr_labels = separate_features_labels(args["train_data"])
    te_features, te_labels = separate_features_labels(args["test_data"])

    classifier = TabPFNClassifier(device='cpu', N_ensemble_configurations=32)
    s_time = time.time()
    classifier.fit(tr_features, tr_labels)
    e_time = time.time()
    print("Time taken by TabPFN for training: {} seconds".format(e_time - s_time))
    y_eval, p_eval = classifier.predict(te_features, return_winning_probability=True)
    print('Accuracy', accuracy_score(te_labels, y_eval))

    te_features["true_labels"] = te_labels
    te_features["pred_labels"] = y_eval
    te_features.to_csv("output_predicted_data.csv", sep="\t", index=None)


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-trdata", "--train_data", required=True, help="Train data")
    arg_parser.add_argument("-tedata", "--test_data", required=True, help="Test data")

    # get argument values
    args = vars(arg_parser.parse_args())
    train_evaluate(args)
