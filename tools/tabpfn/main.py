"""
Tabular data prediction using TabPFN
"""
import argparse
import time

import matplotlib.pyplot as plt
import pandas as pd
from sklearn.metrics import accuracy_score, average_precision_score, precision_recall_curve
from tabpfn import TabPFNClassifier


def separate_features_labels(data):
    df = pd.read_csv(data, sep=",")
    labels = df.iloc[:, -1]
    features = df.iloc[:, :-1]
    return features, labels


def train_evaluate(args):
    """
    Train TabPFN
    """
    tr_features, tr_labels = separate_features_labels(args["train_data"])
    te_features, te_labels = separate_features_labels(args["test_data"])
    classifier = TabPFNClassifier(device='cpu')
    s_time = time.time()
    classifier.fit(tr_features, tr_labels)
    e_time = time.time()
    print("Time taken by TabPFN for training: {} seconds".format(e_time - s_time))
    y_eval = classifier.predict(te_features)
    print('Accuracy', accuracy_score(te_labels, y_eval))
    pred_probas_test = classifier.predict_proba(te_features)
    te_features["predicted_labels"] = y_eval
    te_features.to_csv("output_predicted_data", sep="\t", index=None)
    precision, recall, thresholds = precision_recall_curve(te_labels, pred_probas_test[:, 1])
    average_precision = average_precision_score(te_labels, pred_probas_test[:, 1])
    plt.figure(figsize=(8, 6))
    plt.plot(recall, precision, label=f'Precision-Recall Curve (AP={average_precision:.2f})')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title('Precision-Recall Curve')
    plt.legend(loc='lower left')
    plt.grid(True)
    plt.savefig("output_prec_recall_curve.png")


if __name__ == "__main__":
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument("-trdata", "--train_data", required=True, help="Train data")
    arg_parser.add_argument("-tedata", "--test_data", required=True, help="Test data")
    # get argument values
    args = vars(arg_parser.parse_args())
    train_evaluate(args)
