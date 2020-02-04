import argparse
import json
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import plotly
import plotly.graph_objs as go
import warnings

from keras.models import model_from_json
from keras.utils import plot_model
from sklearn.feature_selection.base import SelectorMixin
from sklearn.metrics import precision_recall_curve, average_precision_score
from sklearn.metrics import roc_curve, auc
from sklearn.pipeline import Pipeline
from galaxy_ml.utils import load_model, read_columns, SafeEval


safe_eval = SafeEval()

# plotly default colors
default_colors = [
    '#1f77b4',  # muted blue
    '#ff7f0e',  # safety orange
    '#2ca02c',  # cooked asparagus green
    '#d62728',  # brick red
    '#9467bd',  # muted purple
    '#8c564b',  # chestnut brown
    '#e377c2',  # raspberry yogurt pink
    '#7f7f7f',  # middle gray
    '#bcbd22',  # curry yellow-green
    '#17becf'   # blue-teal
]


def visualize_pr_curve_plotly(df1, df2, pos_label, title=None):
    """output pr-curve in html using plotly

    df1 : pandas.DataFrame
        Containing y_true
    df2 : pandas.DataFrame
        Containing y_score
    pos_label : None
        The label of positive class
    title : str
        Plot title
    """
    data = []
    for idx in range(df1.shape[1]):
        y_true = df1.iloc[:, idx].values
        y_score = df2.iloc[:, idx].values

        precision, recall, _ = precision_recall_curve(
            y_true, y_score, pos_label=pos_label)
        ap = average_precision_score(
            y_true, y_score, pos_label=pos_label or 1)

        trace = go.Scatter(
            x=recall,
            y=precision,
            mode='lines',
            marker=dict(
                color=default_colors[idx % len(default_colors)]
            ),
            name='%s (area = %.3f)' % (idx, ap)
        )
        data.append(trace)

    layout = go.Layout(
        xaxis=dict(
            title='Recall',
            linecolor='lightslategray',
            linewidth=1
        ),
        yaxis=dict(
            title='Precision',
            linecolor='lightslategray',
            linewidth=1
        ),
        title=dict(
            text=title or 'Precision-Recall Curve',
            x=0.5,
            y=0.92,
            xanchor='center',
            yanchor='top'
        ),
        font=dict(
            family="sans-serif",
            size=11
        ),
        # control backgroud colors
        plot_bgcolor='rgba(255,255,255,0)'
    )
    """
    legend=dict(
        x=0.95,
        y=0,
        traceorder="normal",
        font=dict(
            family="sans-serif",
            size=9,
            color="black"
        ),
        bgcolor="LightSteelBlue",
        bordercolor="Black",
        borderwidth=2
    ),"""

    fig = go.Figure(data=data, layout=layout)

    plotly.offline.plot(fig, filename="output.html", auto_open=False)
    # to be discovered by `from_work_dir`
    os.rename('output.html', 'output')


def visualize_pr_curve_matplotlib(df1, df2, pos_label, title=None):
    """visualize pr-curve using matplotlib and output svg image
    """
    backend = matplotlib.get_backend()
    if "inline" not in backend:
        matplotlib.use("SVG")
    plt.style.use('seaborn-colorblind')
    plt.figure()

    for idx in range(df1.shape[1]):
        y_true = df1.iloc[:, idx].values
        y_score = df2.iloc[:, idx].values

        precision, recall, _ = precision_recall_curve(
            y_true, y_score, pos_label=pos_label)
        ap = average_precision_score(
            y_true, y_score, pos_label=pos_label or 1)

        plt.step(recall, precision, 'r-', color="black", alpha=0.3,
                 lw=1, where="post", label='%s (area = %.3f)' % (idx, ap))

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    title = title or 'Precision-Recall Curve'
    plt.title(title)
    folder = os.getcwd()
    plt.savefig(os.path.join(folder, "output.svg"), format="svg")
    os.rename(os.path.join(folder, "output.svg"),
              os.path.join(folder, "output"))


def visualize_roc_curve_plotly(df1, df2, pos_label,
                               drop_intermediate=True,
                               title=None):
    """output roc-curve in html using plotly

    df1 : pandas.DataFrame
        Containing y_true
    df2 : pandas.DataFrame
        Containing y_score
    pos_label : None
        The label of positive class
    drop_intermediate : bool
        Whether to drop some suboptimal thresholds
    title : str
        Plot title
    """
    data = []
    for idx in range(df1.shape[1]):
        y_true = df1.iloc[:, idx].values
        y_score = df2.iloc[:, idx].values

        fpr, tpr, _ = roc_curve(y_true, y_score, pos_label=pos_label,
                                drop_intermediate=drop_intermediate)
        roc_auc = auc(fpr, tpr)

        trace = go.Scatter(
            x=fpr,
            y=tpr,
            mode='lines',
            marker=dict(
                color=default_colors[idx % len(default_colors)]
            ),
            name='%s (area = %.3f)' % (idx, roc_auc)
        )
        data.append(trace)

    layout = go.Layout(
        xaxis=dict(
            title='False Positive Rate',
            linecolor='lightslategray',
            linewidth=1
        ),
        yaxis=dict(
            title='True Positive Rate',
            linecolor='lightslategray',
            linewidth=1
        ),
        title=dict(
            text=title or 'Receiver Operating Characteristic (ROC) Curve',
            x=0.5,
            y=0.92,
            xanchor='center',
            yanchor='top'
        ),
        font=dict(
            family="sans-serif",
            size=11
        ),
        # control backgroud colors
        plot_bgcolor='rgba(255,255,255,0)'
    )
    """
    # legend=dict(
            # x=0.95,
            # y=0,
            # traceorder="normal",
            # font=dict(
            #    family="sans-serif",
            #    size=9,
            #    color="black"
            # ),
            # bgcolor="LightSteelBlue",
            # bordercolor="Black",
            # borderwidth=2
        # ),
    """

    fig = go.Figure(data=data, layout=layout)

    plotly.offline.plot(fig, filename="output.html", auto_open=False)
    # to be discovered by `from_work_dir`
    os.rename('output.html', 'output')


def visualize_roc_curve_matplotlib(df1, df2, pos_label,
                                   drop_intermediate=True,
                                   title=None):
    """visualize roc-curve using matplotlib and output svg image
    """
    backend = matplotlib.get_backend()
    if "inline" not in backend:
        matplotlib.use("SVG")
    plt.style.use('seaborn-colorblind')
    plt.figure()

    for idx in range(df1.shape[1]):
        y_true = df1.iloc[:, idx].values
        y_score = df2.iloc[:, idx].values

        fpr, tpr, _ = roc_curve(y_true, y_score, pos_label=pos_label,
                                drop_intermediate=drop_intermediate)
        roc_auc = auc(fpr, tpr)

        plt.step(fpr, tpr, 'r-', color="black", alpha=0.3, lw=1,
                 where="post", label='%s (area = %.3f)' % (idx, roc_auc))

    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    title = title or 'Receiver Operating Characteristic (ROC) Curve'
    plt.title(title)
    folder = os.getcwd()
    plt.savefig(os.path.join(folder, "output.svg"), format="svg")
    os.rename(os.path.join(folder, "output.svg"),
              os.path.join(folder, "output"))


def main(inputs, infile_estimator=None, infile1=None,
         infile2=None, outfile_result=None,
         outfile_object=None, groups=None,
         ref_seq=None, intervals=None,
         targets=None, fasta_path=None,
         model_config=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_estimator : str, default is None
        File path to estimator

    infile1 : str, default is None
        File path to dataset containing features or true labels.

    infile2 : str, default is None
        File path to dataset containing target values or predicted
        probabilities.

    outfile_result : str, default is None
        File path to save the results, either cv_results or test result

    outfile_object : str, default is None
        File path to save searchCV object

    groups : str, default is None
        File path to dataset containing groups labels

    ref_seq : str, default is None
        File path to dataset containing genome sequence file

    intervals : str, default is None
        File path to dataset containing interval file

    targets : str, default is None
        File path to dataset compressed target bed file

    fasta_path : str, default is None
        File path to dataset containing fasta file

    model_config : str, default is None
        File path to dataset containing JSON config for neural networks
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    title = params['plotting_selection']['title'].strip()
    plot_type = params['plotting_selection']['plot_type']
    plot_format = params['plotting_selection']['plot_format']

    if plot_type == 'feature_importances':
        with open(infile_estimator, 'rb') as estimator_handler:
            estimator = load_model(estimator_handler)

        column_option = (params['plotting_selection']
                               ['column_selector_options']
                               ['selected_column_selector_option'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = (params['plotting_selection']
                       ['column_selector_options']['col1'])
        else:
            c = None

        _, input_df = read_columns(infile1, c=c,
                                   c_option=column_option,
                                   return_df=True,
                                   sep='\t', header='infer',
                                   parse_dates=True)

        feature_names = input_df.columns.values

        if isinstance(estimator, Pipeline):
            for st in estimator.steps[:-1]:
                if isinstance(st[-1], SelectorMixin):
                    mask = st[-1].get_support()
                    feature_names = feature_names[mask]
            estimator = estimator.steps[-1][-1]

        if hasattr(estimator, 'coef_'):
            coefs = estimator.coef_
        else:
            coefs = getattr(estimator, 'feature_importances_', None)
        if coefs is None:
            raise RuntimeError('The classifier does not expose '
                               '"coef_" or "feature_importances_" '
                               'attributes')

        threshold = params['plotting_selection']['threshold']
        if threshold is not None:
            mask = (coefs > threshold) | (coefs < -threshold)
            coefs = coefs[mask]
            feature_names = feature_names[mask]

        # sort
        indices = np.argsort(coefs)[::-1]

        trace = go.Bar(x=feature_names[indices],
                       y=coefs[indices])
        layout = go.Layout(title=title or "Feature Importances")
        fig = go.Figure(data=[trace], layout=layout)

        plotly.offline.plot(fig, filename="output.html",
                            auto_open=False)
        # to be discovered by `from_work_dir`
        os.rename('output.html', 'output')

        return 0

    elif plot_type in ('pr_curve', 'roc_curve'):
        df1 = pd.read_csv(infile1, sep='\t', header='infer')
        df2 = pd.read_csv(infile2, sep='\t', header='infer').astype(np.float32)

        minimum = params['plotting_selection']['report_minimum_n_positives']
        # filter out columns whose n_positives is beblow the threhold
        if minimum:
            mask = df1.sum(axis=0) >= minimum
            df1 = df1.loc[:, mask]
            df2 = df2.loc[:, mask]

        pos_label = params['plotting_selection']['pos_label'].strip() \
            or None

        if plot_type == 'pr_curve':
            if plot_format == 'plotly_html':
                visualize_pr_curve_plotly(df1, df2, pos_label, title=title)
            else:
                visualize_pr_curve_matplotlib(df1, df2, pos_label, title)
        else:          # 'roc_curve'
            drop_intermediate = (params['plotting_selection']
                                       ['drop_intermediate'])
            if plot_format == 'plotly_html':
                visualize_roc_curve_plotly(df1, df2, pos_label,
                                           drop_intermediate=drop_intermediate,
                                           title=title)
            else:
                visualize_roc_curve_matplotlib(
                    df1, df2, pos_label,
                    drop_intermediate=drop_intermediate,
                    title=title)

        return 0

    elif plot_type == 'rfecv_gridscores':
        input_df = pd.read_csv(infile1, sep='\t', header='infer')
        scores = input_df.iloc[:, 0]
        steps = params['plotting_selection']['steps'].strip()
        steps = safe_eval(steps)

        data = go.Scatter(
            x=list(range(len(scores))),
            y=scores,
            text=[str(_) for _ in steps] if steps else None,
            mode='lines'
        )
        layout = go.Layout(
            xaxis=dict(title="Number of features selected"),
            yaxis=dict(title="Cross validation score"),
            title=dict(
                text=title or None,
                x=0.5,
                y=0.92,
                xanchor='center',
                yanchor='top'
            ),
            font=dict(
                family="sans-serif",
                size=11
            ),
            # control backgroud colors
            plot_bgcolor='rgba(255,255,255,0)'
        )
        """
        # legend=dict(
                # x=0.95,
                # y=0,
                # traceorder="normal",
                # font=dict(
                #    family="sans-serif",
                #    size=9,
                #    color="black"
                # ),
                # bgcolor="LightSteelBlue",
                # bordercolor="Black",
                # borderwidth=2
            # ),
        """

        fig = go.Figure(data=[data], layout=layout)
        plotly.offline.plot(fig, filename="output.html",
                            auto_open=False)
        # to be discovered by `from_work_dir`
        os.rename('output.html', 'output')

        return 0

    elif plot_type == 'learning_curve':
        input_df = pd.read_csv(infile1, sep='\t', header='infer')
        plot_std_err = params['plotting_selection']['plot_std_err']
        data1 = go.Scatter(
            x=input_df['train_sizes_abs'],
            y=input_df['mean_train_scores'],
            error_y=dict(
                array=input_df['std_train_scores']
            ) if plot_std_err else None,
            mode='lines',
            name="Train Scores",
        )
        data2 = go.Scatter(
            x=input_df['train_sizes_abs'],
            y=input_df['mean_test_scores'],
            error_y=dict(
                array=input_df['std_test_scores']
            ) if plot_std_err else None,
            mode='lines',
            name="Test Scores",
        )
        layout = dict(
            xaxis=dict(
                title='No. of samples'
            ),
            yaxis=dict(
                title='Performance Score'
            ),
            # modify these configurations to customize image
            title=dict(
                text=title or 'Learning Curve',
                x=0.5,
                y=0.92,
                xanchor='center',
                yanchor='top'
            ),
            font=dict(
                family="sans-serif",
                size=11
            ),
            # control backgroud colors
            plot_bgcolor='rgba(255,255,255,0)'
        )
        """
        # legend=dict(
                # x=0.95,
                # y=0,
                # traceorder="normal",
                # font=dict(
                #    family="sans-serif",
                #    size=9,
                #    color="black"
                # ),
                # bgcolor="LightSteelBlue",
                # bordercolor="Black",
                # borderwidth=2
            # ),
        """

        fig = go.Figure(data=[data1, data2], layout=layout)
        plotly.offline.plot(fig, filename="output.html",
                            auto_open=False)
        # to be discovered by `from_work_dir`
        os.rename('output.html', 'output')

        return 0

    elif plot_type == 'keras_plot_model':
        with open(model_config, 'r') as f:
            model_str = f.read()
        model = model_from_json(model_str)
        plot_model(model, to_file="output.png")
        os.rename('output.png', 'output')

        return 0

    # save pdf file to disk
    # fig.write_image("image.pdf", format='pdf')
    # fig.write_image("image.pdf", format='pdf', width=340*2, height=226*2)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--estimator", dest="infile_estimator")
    aparser.add_argument("-X", "--infile1", dest="infile1")
    aparser.add_argument("-y", "--infile2", dest="infile2")
    aparser.add_argument("-O", "--outfile_result", dest="outfile_result")
    aparser.add_argument("-o", "--outfile_object", dest="outfile_object")
    aparser.add_argument("-g", "--groups", dest="groups")
    aparser.add_argument("-r", "--ref_seq", dest="ref_seq")
    aparser.add_argument("-b", "--intervals", dest="intervals")
    aparser.add_argument("-t", "--targets", dest="targets")
    aparser.add_argument("-f", "--fasta_path", dest="fasta_path")
    aparser.add_argument("-c", "--model_config", dest="model_config")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile1, args.infile2,
         args.outfile_result, outfile_object=args.outfile_object,
         groups=args.groups, ref_seq=args.ref_seq, intervals=args.intervals,
         targets=args.targets, fasta_path=args.fasta_path,
         model_config=args.model_config)
