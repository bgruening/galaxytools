import argparse
import pandas as pd
import numpy as np
import plotly
import plotly.graph_objs as go


def main(infile_input, infile_output):
    """
    Produce an interactive actual vs predicted curves and residual plots
    Args:
        infile_input: str, input tabular file with true values
        infile_output: str, input tabular file with predicted values
    """

    df_input = pd.read_csv(infile_input, sep='\t', parse_dates=True)
    df_output = pd.read_csv(infile_output, sep='\t', parse_dates=True)
    true_values = df_input.iloc[:, -1].copy()
    predicted_values = df_output.iloc[:, -1].copy()
    axis_labels = list(range(1, len(true_values)+1))

    # true vs predicted curves
    trace_true = go.Scatter(
        x=axis_labels,
        y=true_values,
        mode='lines+markers',
        name='True values'
    )

    trace_predicted = go.Scatter(
        x=axis_labels,
        y=predicted_values,
        mode='lines+markers',
        name='Predicted values'
    )

    layout_tp = go.Layout(
        title='True vs predicted values',
        xaxis=dict(title='Number of data points'),
        yaxis=dict(title='Values')
    )

    data_tp = [trace_true, trace_predicted]
    fig_tp = go.Figure(data=data_tp, layout=layout_tp)
    plotly.offline.plot(fig_tp, filename="output_actual_vs_pred.html", auto_open=False)

    # scatter plot
    max_tv = int(max(true_values))
    x_y_values = list(range(0, max_tv))

    true_mean = np.mean(true_values)
    res_true_predicted = np.sum((true_values - predicted_values) ** 2)
    res_total = np.sum((true_values - true_mean) ** 2)
    r2 = 1 - (res_true_predicted / float(res_total))
    rmse = np.sqrt(np.mean([(x - y) ** 2 for x, y in zip(true_values, predicted_values)]))

    trace_x_eq_y = go.Scatter(
        x=x_y_values,
        y=x_y_values,
        mode='lines',
        name='X = Y curve'
    )

    trace_true_pred = go.Scatter(
        x=true_values,
        y=predicted_values,
        mode='markers',
        name='True and predicted values'
    )

    layout_true_pred = go.Layout(
        title='True vs predicted values (RMSE: %s, R2: %s)' % (str(np.round(rmse, 2)), str(np.round(r2, 2))),
        xaxis=dict(title='True values'),
        yaxis=dict(title='Predicted values')
    )

    data_true_pred = [trace_true_pred, trace_x_eq_y]
    fig_true_pred = go.Figure(data=data_true_pred, layout=layout_true_pred)
    plotly.offline.plot(fig_true_pred, filename="output_scatter_plot.html", auto_open=False)

    # residual plot
    residual = predicted_values - true_values
    trace_residual = go.Scatter(
        x=predicted_values,
        y=residual,
        mode='markers'
    )

    layout_residual = go.Layout(
        title='Residual vs predicted values',
        xaxis=dict(title='Predicted values'),
        yaxis=dict(title='Residual (Predicted - True)')
    )

    data_residual = [trace_residual]
    fig_residual = go.Figure(data=data_residual, layout=layout_residual)
    plotly.offline.plot(fig_residual, filename="output_residual_plot.html", auto_open=False)


if __name__ == "__main__":
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--input", dest="infile_input", required=True)
    aparser.add_argument("-j", "--output", dest="infile_output", required=True)
    args = aparser.parse_args()
    main(args.infile_input, args.infile_output)
