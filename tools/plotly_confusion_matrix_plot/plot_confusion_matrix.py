import argparse
import pandas as pd
import plotly
import pickle
import plotly.graph_objs as go
from sklearn.metrics import confusion_matrix, precision_recall_fscore_support, roc_curve, auc
from scipy import interp
from sklearn.preprocessing import label_binarize
from sklearn.model_selection import StratifiedKFold


def main(infile_input, infile_output, infile_trained_model):
    """
    Produce an interactive confusion matrix (heatmap) plot
    Args:
        infile_input: str, input tabular file
        infile_output: str, output tabular file
    """

    df_input = pd.read_csv(infile_input, sep='\t', parse_dates=True)
    df_output = pd.read_csv(infile_output, sep='\t', parse_dates=True)
    true_labels = df_input.iloc[:, -1].copy()
    predicted_labels = df_output.iloc[:, -1].copy()
    axis_labels = list(set(true_labels))
    c_matrix = confusion_matrix(true_labels, predicted_labels)
    data = [
        go.Heatmap(
            z=c_matrix,
            x=axis_labels,
            y=axis_labels,
            colorscale='Portland',
        )
    ]

    layout = go.Layout(
        title='Confusion Matrix between true and predicted class labels',
        xaxis=dict(title='Predicted class labels'),
        yaxis=dict(title='True class labels')
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename="output_confusion.html", auto_open=False)
    
    # plot precision, recall and f_score for each class label
    precision, recall, f_score, _ = precision_recall_fscore_support(true_labels, predicted_labels)
    
    trace_precision = go.Scatter(
        x = axis_labels,
        y = precision,
        mode = 'lines+markers',
        name = 'Precision'
    )
    
    trace_recall = go.Scatter(
        x = axis_labels,
        y = recall,
        mode = 'lines+markers',
        name = 'Recall'
    )
    
    trace_fscore = go.Scatter(
        x = axis_labels,
        y = f_score,
        mode = 'lines+markers',
        name = 'F-score'
    )
    
    layout_prf = go.Layout(
        title='Precision, recall and f-score of true and predicted class labels',
        xaxis=dict(title='Class labels'),
        yaxis=dict(title='Precision, recall and f-score')
    )
    
    data_prf = [trace_precision, trace_recall, trace_fscore]
    fig_prf = go.Figure(data=data_prf, layout=layout_prf)
    plotly.offline.plot(fig_prf, filename="output_prf.html", auto_open=False)
    
    # plot roc and auc curves for different classes
    with open(infile_trained_model, 'rb') as model_file:
        model = pickle.load(model_file)
    test_data = df_input.drop("target", axis=1)
    y_score = model.predict_proba(test_data)
    fpr = dict()
    tpr = dict()
    roc_auc = dict()
    true_labels_list = true_labels.tolist()
    one_hot_labels = label_binarize(true_labels_list, classes=axis_labels)
    for i in axis_labels:
        fpr[i], tpr[i], _ = roc_curve(one_hot_labels[:, i], y_score[:, i])
        roc_auc[i] = auc(fpr[i], tpr[i])
    fpr["micro"], tpr["micro"], _ = roc_curve(one_hot_labels.ravel(), y_score.ravel())
    roc_auc["micro"] = auc(fpr["micro"], tpr["micro"])
    
    data_roc = list()
    for i in range(len(axis_labels)):
        trace = go.Scatter(
            x = fpr[i],
            y = tpr[i],
            mode = 'lines+markers',
            name = 'ROC curve of class {0} (AUC = {1:0.2f})'.format(i, roc_auc[i])
        )
        data_roc.append(trace)
    
    trace_diag = go.Scatter(
        x = [0, 1],
        y = [0, 1],
        mode = 'lines',
        name = 'Chance'
    )
    
    layout_roc = go.Layout(
        title='Receiver operating characteristics (ROC) and area under curve (AUC)',
        xaxis=dict(title='False positive rate'),
        yaxis=dict(title='True positive rate')
    )
    data_roc.append(trace_diag)
    fig_roc = go.Figure(data=data_roc, layout=layout_roc)
    plotly.offline.plot(fig_roc, filename="output_roc.html", auto_open=False)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--input", dest="infile_input", required=True)
    aparser.add_argument("-j", "--output", dest="infile_output", required=True)
    aparser.add_argument("-k", "--model", dest="infile_trained_model", required=True)
    args = aparser.parse_args()
    main(args.infile_input, args.infile_output, args.infile_trained_model)
