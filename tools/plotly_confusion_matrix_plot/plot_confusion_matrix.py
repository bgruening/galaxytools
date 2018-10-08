import sys
import argparse
import pandas as pd
import plotly
import plotly.graph_objs as go
from sklearn.metrics import confusion_matrix


def main(infile_input, infile_output):
    """
    Produce an interactive confusion matrix (heatmap) plot
    Args:
        infile_input: str, input tabular file
        infile_output: str, output tabular file
    """
    
    df_input = pd.read_csv(infile_input, sep='\t', parse_dates=True)
    df_output = pd.read_csv(infile_output, sep='\t', parse_dates=True)
    true_labels = df_input.iloc[:,-1].copy()
    predicted_labels = df_output.iloc[:,-1].copy()
    
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
        xaxis = dict(title='Predicted class labels'),
        yaxis = dict(title='True class labels')
    )

    fig = go.Figure(data=data, layout=layout)
    plotly.offline.plot(fig, filename="output", auto_open=False)

if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument( "-i", "--input", dest="infile_input", required=True)
    aparser.add_argument( "-j", "--output", dest="infile_output", required=True)
    args = aparser.parse_args()

    main(args.infile_input, args.infile_output)
