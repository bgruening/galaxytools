import argparse
import json
import pandas as pd
import warnings

from scipy.io import mmread
from sklearn.pipeline import Pipeline
from sklearn.metrics.scorer import _check_multimetric_scoring
from sklearn.model_selection._validation import _score
from galaxy_ml.utils import get_scoring, load_model, read_columns


N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))


def _get_X_y(params, infile1, infile2):
    """ read from inputs and output X and y

    Parameters
    ----------
    params : dict
        Tool inputs parameter
    infile1 : str
        File path to dataset containing features
    infile2 : str
        File path to dataset containing target values

    """
    # store read dataframe object
    loaded_df = {}

    input_type = params['input_options']['selected_input']
    # tabular input
    if input_type == 'tabular':
        header = 'infer' if params['input_options']['header1'] else None
        column_option = (params['input_options']['column_selector_options_1']
                         ['selected_column_selector_option'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = params['input_options']['column_selector_options_1']['col1']
        else:
            c = None

        df_key = infile1 + repr(header)
        df = pd.read_csv(infile1, sep='\t', header=header,
                         parse_dates=True)
        loaded_df[df_key] = df

        X = read_columns(df, c=c, c_option=column_option).astype(float)
    # sparse input
    elif input_type == 'sparse':
        X = mmread(open(infile1, 'r'))

    # Get target y
    header = 'infer' if params['input_options']['header2'] else None
    column_option = (params['input_options']['column_selector_options_2']
                     ['selected_column_selector_option2'])
    if column_option in ['by_index_number', 'all_but_by_index_number',
                         'by_header_name', 'all_but_by_header_name']:
        c = params['input_options']['column_selector_options_2']['col2']
    else:
        c = None

    df_key = infile2 + repr(header)
    if df_key in loaded_df:
        infile2 = loaded_df[df_key]
    else:
        infile2 = pd.read_csv(infile2, sep='\t',
                              header=header, parse_dates=True)
        loaded_df[df_key] = infile2

    y = read_columns(
            infile2,
            c=c,
            c_option=column_option,
            sep='\t',
            header=header,
            parse_dates=True)
    if len(y.shape) == 2 and y.shape[1] == 1:
        y = y.ravel()

    return X, y


def main(inputs, infile_estimator, outfile_eval,
         infile_weights=None, infile1=None,
         infile2=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_estimator : strgit
        File path to trained estimator input

    outfile_eval : str
        File path to save the evalulation results, tabular

    infile_weights : str
        File path to weights input

    infile1 : str
        File path to dataset containing features

    infile2 : str
        File path to dataset containing target values
    """
    warnings.filterwarnings('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    X_test, y_test = _get_X_y(params, infile1, infile2)

    # load model
    with open(infile_estimator, 'rb') as est_handler:
        estimator = load_model(est_handler)

    main_est = estimator
    if isinstance(estimator, Pipeline):
        main_est = estimator.steps[-1][-1]
    if hasattr(main_est, 'config') and hasattr(main_est, 'load_weights'):
        if not infile_weights or infile_weights == 'None':
            raise ValueError("The selected model skeleton asks for weights, "
                             "but dataset for weights wan not selected!")
        main_est.load_weights(infile_weights)

    # handle scorer, convert to scorer dict
    scoring = params['scoring']
    scorer = get_scoring(scoring)
    scorer, _ = _check_multimetric_scoring(estimator, scoring=scorer)

    if hasattr(estimator, 'evaluate'):
        scores = estimator.evaluate(X_test, y_test=y_test,
                                    scorer=scorer,
                                    is_multimetric=True)
    else:
        scores = _score(estimator, X_test, y_test, scorer,
                        is_multimetric=True)

    # handle output
    for name, score in scores.items():
        scores[name] = [score]
    df = pd.DataFrame(scores)
    df = df[sorted(df.columns)]
    df.to_csv(path_or_buf=outfile_eval, sep='\t',
              header=True, index=False)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--infile_estimator", dest="infile_estimator")
    aparser.add_argument("-w", "--infile_weights", dest="infile_weights")
    aparser.add_argument("-X", "--infile1", dest="infile1")
    aparser.add_argument("-y", "--infile2", dest="infile2")
    aparser.add_argument("-O", "--outfile_eval", dest="outfile_eval")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.outfile_eval,
         infile_weights=args.infile_weights, infile1=args.infile1,
         infile2=args.infile2)
