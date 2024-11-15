import argparse

import matplotlib.pyplot as plt
import pandas as pd
import plotly
import plotly.graph_objs as go
from galaxy_ml.model_persist import load_model_from_h5
from galaxy_ml.utils import clean_params
from sklearn.metrics import (
    auc,
    confusion_matrix,
    precision_recall_fscore_support,
    roc_curve,
)
from sklearn.preprocessing import label_binarize


def main(infile_input, infile_output, infile_trained_model):
    """
    Produce an interactive confusion matrix (heatmap), precision, recall, fscore and auc plots
    Args:
        infile_input: str, input tabular file with true labels
        infile_output: str, input tabular file with predicted labels
        infile_trained_model: str, input trained model file (h5mlm)
    """

    df_input = pd.read_csv(infile_input, sep="\t", parse_dates=True)
    df_output = pd.read_csv(infile_output, sep="\t", parse_dates=True)
    true_labels = df_input.iloc[:, -1].copy()
    predicted_labels = df_output.iloc[:, -1].copy()
    axis_labels = list(set(true_labels))
    c_matrix = confusion_matrix(true_labels, predicted_labels)
    fig, ax = plt.subplots(figsize=(7, 7))
    im = plt.imshow(c_matrix, cmap="viridis")
    # add number of samples to each cell of confusion matrix plot
    for i in range(len(c_matrix)):
        for j in range(len(c_matrix)):
            ax.text(j, i, c_matrix[i, j], ha="center", va="center", color="k")
    ax.set_ylabel("True class labels")
    ax.set_xlabel("Predicted class labels")
    ax.set_title("Confusion Matrix between true and predicted class labels")
    ax.set_xticks(axis_labels)
    ax.set_yticks(axis_labels)
    fig.colorbar(im, ax=ax)
    fig.tight_layout()
    plt.savefig("output_confusion.png", dpi=120)

    # plot precision, recall and f_score for each class label
    precision, recall, f_score, _ = precision_recall_fscore_support(
        true_labels, predicted_labels
    )

    trace_precision = go.Scatter(
        x=axis_labels, y=precision, mode="lines+markers", name="Precision"
    )

    trace_recall = go.Scatter(
        x=axis_labels, y=recall, mode="lines+markers", name="Recall"
    )

    trace_fscore = go.Scatter(
        x=axis_labels, y=f_score, mode="lines+markers", name="F-score"
    )

    layout_prf = go.Layout(
        title="Precision, recall and f-score of true and predicted class labels",
        xaxis=dict(title="Class labels"),
        yaxis=dict(title="Precision, recall and f-score"),
    )

    data_prf = [trace_precision, trace_recall, trace_fscore]
    fig_prf = go.Figure(data=data_prf, layout=layout_prf)
    plotly.offline.plot(fig_prf, filename="output_prf.html", auto_open=False)

    # plot roc and auc curves for different classes
    classifier_object = load_model_from_h5(infile_trained_model)
    model = clean_params(classifier_object)

    # remove the last column (label column)
    test_data = df_input.iloc[:, :-1]
    model_items = dir(model)

    try:
        # find the probability estimating method
        if "predict_proba" in model_items:
            y_score = model.predict_proba(test_data)
        elif "decision_function" in model_items:
            y_score = model.decision_function(test_data)

        true_labels_list = true_labels.tolist()
        one_hot_labels = label_binarize(true_labels_list, classes=axis_labels)
        data_roc = list()

        if len(axis_labels) > 2:
            fpr = dict()
            tpr = dict()
            roc_auc = dict()
            for i in axis_labels:
                fpr[i], tpr[i], _ = roc_curve(one_hot_labels[:, i], y_score[:, i])
                roc_auc[i] = auc(fpr[i], tpr[i])
            for i in range(len(axis_labels)):
                trace = go.Scatter(
                    x=fpr[i],
                    y=tpr[i],
                    mode="lines+markers",
                    name="ROC curve of class {0} (AUC = {1:0.2f})".format(
                        i, roc_auc[i]
                    ),
                )
                data_roc.append(trace)
        else:
            try:
                y_score_binary = y_score[:, 1]
            except Exception:
                y_score_binary = y_score
            fpr, tpr, _ = roc_curve(one_hot_labels, y_score_binary, pos_label=1)
            roc_auc = auc(fpr, tpr)
            trace = go.Scatter(
                x=fpr,
                y=tpr,
                mode="lines+markers",
                name="ROC curve (AUC = {0:0.2f})".format(roc_auc),
            )
            data_roc.append(trace)

        trace_diag = go.Scatter(x=[0, 1], y=[0, 1], mode="lines", name="Chance")
        data_roc.append(trace_diag)
        layout_roc = go.Layout(
            title="Receiver operating characteristics (ROC) and area under curve (AUC)",
            xaxis=dict(title="False positive rate"),
            yaxis=dict(title="True positive rate"),
        )

        fig_roc = go.Figure(data=data_roc, layout=layout_roc)
        plotly.offline.plot(fig_roc, filename="output_roc.html", auto_open=False)

    except Exception as exp:
        print(
            "Plotting the ROC-AUC graph failed. This exception was raised: {}".format(
                exp
            )
        )


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--input", dest="infile_input", required=True)
    aparser.add_argument("-j", "--output", dest="infile_output", required=True)
    aparser.add_argument("-k", "--model", dest="infile_trained_model", required=True)
    args = aparser.parse_args()
    main(args.infile_input, args.infile_output, args.infile_trained_model)
