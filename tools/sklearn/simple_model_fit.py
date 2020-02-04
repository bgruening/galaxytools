import argparse
import json
import pandas as pd
import pickle

from galaxy_ml.utils import load_model, read_columns
from sklearn.pipeline import Pipeline


N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))


# TODO import from galaxy_ml.utils in future versions
def clean_params(estimator, n_jobs=None):
    """clean unwanted hyperparameter settings

    If n_jobs is not None, set it into the estimator, if applicable

    Return
    ------
    Cleaned estimator object
    """
    ALLOWED_CALLBACKS = ('EarlyStopping', 'TerminateOnNaN',
                         'ReduceLROnPlateau', 'CSVLogger', 'None')

    estimator_params = estimator.get_params()

    for name, p in estimator_params.items():
        # all potential unauthorized file write
        if name == 'memory' or name.endswith('__memory') \
                or name.endswith('_path'):
            new_p = {name: None}
            estimator.set_params(**new_p)
        elif n_jobs is not None and (name == 'n_jobs' or
                                     name.endswith('__n_jobs')):
            new_p = {name: n_jobs}
            estimator.set_params(**new_p)
        elif name.endswith('callbacks'):
            for cb in p:
                cb_type = cb['callback_selection']['callback_type']
                if cb_type not in ALLOWED_CALLBACKS:
                    raise ValueError(
                        "Prohibited callback type: %s!" % cb_type)

    return estimator


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


def main(inputs, infile_estimator, infile1, infile2, out_object,
         out_weights=None):
    """ main

    Parameters
    ----------
    inputs : str
        File path to galaxy tool parameter

    infile_estimator : str
        File paths of input estimator

    infile1 : str
        File path to dataset containing features

    infile2 : str
        File path to dataset containing target labels

    out_object : str
        File path for output of fitted model or skeleton

    out_weights : str
        File path for output of weights

    """
    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    # load model
    with open(infile_estimator, 'rb') as est_handler:
        estimator = load_model(est_handler)
    estimator = clean_params(estimator, n_jobs=N_JOBS)

    X_train, y_train = _get_X_y(params, infile1, infile2)

    estimator.fit(X_train, y_train)
    
    main_est = estimator
    if isinstance(main_est, Pipeline):
        main_est = main_est.steps[-1][-1]
    if hasattr(main_est, 'model_') \
            and hasattr(main_est, 'save_weights'):
        if out_weights:
            main_est.save_weights(out_weights)
        del main_est.model_
        del main_est.fit_params
        del main_est.model_class_
        del main_est.validation_data
        if getattr(main_est, 'data_generator_', None):
            del main_est.data_generator_

    with open(out_object, 'wb') as output_handler:
        pickle.dump(estimator, output_handler,
                    pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-X", "--infile_estimator", dest="infile_estimator")
    aparser.add_argument("-y", "--infile1", dest="infile1")
    aparser.add_argument("-g", "--infile2", dest="infile2")
    aparser.add_argument("-o", "--out_object", dest="out_object")
    aparser.add_argument("-t", "--out_weights", dest="out_weights")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile1,
         args.infile2, args.out_object, args.out_weights)
