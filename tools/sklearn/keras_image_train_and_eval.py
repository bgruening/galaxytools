import argparse
import json
import numpy as np
import pandas as pd
import pickle
import warnings

from itertools import chain
from sklearn.metrics.scorer import _check_multimetric_scoring
from sklearn.utils import indexable, safe_indexing
from galaxy_ml.model_validations import train_test_split
from galaxy_ml.keras_galaxy_models import (_predict_generator,
                                           KerasGBatchClassifier)
from galaxy_ml.preprocessors import ImageDataFrameBatchGenerator
from galaxy_ml.utils import (SafeEval, get_scoring, load_model,
                             read_columns, clean_params,
                             get_main_estimator, gen_compute_scores)


WORKING_DIR = __import__('os').getcwd()
IMAGES_DIR = __import__('os').path.join(WORKING_DIR, 'images')
NON_SEARCHABLE = ('n_jobs', 'pre_dispatch', 'memory', '_path', '_dir',
                  'nthread', 'callbacks')
ALLOWED_CALLBACKS = ('EarlyStopping', 'TerminateOnNaN', 'ReduceLROnPlateau',
                     'CSVLogger', 'None')


def _eval_swap_params(params_builder):
    swap_params = {}

    for p in params_builder['param_set']:
        swap_value = p['sp_value'].strip()
        if swap_value == '':
            continue

        param_name = p['sp_name']
        if param_name.lower().endswith(NON_SEARCHABLE):
            warnings.warn("Warning: `%s` is not eligible for search and was "
                          "omitted!" % param_name)
            continue

        if not swap_value.startswith(':'):
            safe_eval = SafeEval(load_scipy=True, load_numpy=True)
            ev = safe_eval(swap_value)
        else:
            # Have `:` before search list, asks for estimator evaluatio
            safe_eval_es = SafeEval(load_estimators=True)
            swap_value = swap_value[1:].strip()
            # TODO maybe add regular express check
            ev = safe_eval_es(swap_value)

        swap_params[param_name] = ev

    return swap_params


def train_test_split_none(*arrays, **kwargs):
    """extend train_test_split to take None arrays
    and support split by group names.
    """
    nones = []
    new_arrays = []
    for idx, arr in enumerate(arrays):
        if arr is None:
            nones.append(idx)
        else:
            new_arrays.append(arr)

    if kwargs['shuffle'] == 'None':
        kwargs['shuffle'] = None

    group_names = kwargs.pop('group_names', None)

    if group_names is not None and group_names.strip():
        group_names = [name.strip() for name in
                       group_names.split(',')]
        new_arrays = indexable(*new_arrays)
        groups = kwargs['labels']
        n_samples = new_arrays[0].shape[0]
        index_arr = np.arange(n_samples)
        test = index_arr[np.isin(groups, group_names)]
        train = index_arr[~np.isin(groups, group_names)]
        rval = list(chain.from_iterable(
            (safe_indexing(a, train),
             safe_indexing(a, test)) for a in new_arrays))
    else:
        rval = train_test_split(*new_arrays, **kwargs)

    for pos in nones:
        rval[pos * 2: 2] = [None, None]

    return rval


def _handle_image_generator_params(params, image_df):
    """reconstruct generator kwargs from tool inputs
    """
    safe_eval = SafeEval()

    options = {}

    headers = image_df.columns
    options['x_col'] = headers[params['x_col'][0] - 1]

    y_col = list(map(lambda x: x-1, params['y_col']))
    if len(y_col) == 1:
        options['y_col'] = headers[y_col[0]]
    else:
        options['y_col'] = list(headers[y_col])

    weight_col = params['weight_col'][0]
    if weight_col is None:
        options['weight_col'] = None
    else:
        options['weight_col'] = headers[weight_col - 1]

    other_options = params['options']
    for k, v in other_options.items():
        if k == 'target_size' or k.endswith('_range'):
            v = v.strip()
            if not v:
                other_options[k] = None
            else:
                other_options[k] = safe_eval(v)
        if k == 'classes':
            v = v.strip()
            if not v:
                other_options[k] = None
            else:
                other_options[k] = [x.strip() for x in v.split(',')]

    options.update(other_options)

    return options


def main(inputs, infile_estimator, infile_images, infile_dataframe,
         outfile_result, outfile_object=None, outfile_weights=None,
         outfile_y_true=None, outfile_y_preds=None, groups=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_estimator : str
        File path to estimator

    infile_images : str
        File path to datasets containing images

    infile_dataframe : str
        File path to tabular dataset containing image information

    outfile_result : str
        File path to save the results, either cv_results or test result

    outfile_object : str, optional
        File path to save searchCV object

    outfile_weights : str, optional
        File path to save deep learning model weights

    outfile_y_true : str, optional
        File path to target values for prediction

    outfile_y_preds : str, optional
        File path to save deep learning model weights

    groups : str
        File path to dataset containing groups labels
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    #  load estimator
    with open(infile_estimator, 'rb') as estimator_handler:
        estimator = load_model(estimator_handler)

    estimator = clean_params(estimator)

    if not isinstance(estimator, KerasGBatchClassifier):
        raise ValueError(
            "Only `galaxy_ml.keras_galaxy_models.KerasGBatchClassifier` "
            "is supported!")

    # swap hyperparameter
    swapping = params['experiment_schemes']['hyperparams_swapping']
    swap_params = _eval_swap_params(swapping)
    estimator.set_params(**swap_params)

    # read DataFrame for images
    data_frame = pd.read_csv(infile_dataframe, sep='\t', header='infer')

    kwargs = _handle_image_generator_params(params['input_options'],
                                            data_frame)

    # build data generator
    image_generator = ImageDataFrameBatchGenerator(dataframe=data_frame,
                                                   directory=IMAGES_DIR,
                                                   **kwargs)
    estimator.set_params(data_batch_generator=image_generator)

    # Get X and y
    X = np.arange(data_frame.shape[0])[:, np.newaxis]

    if isinstance(kwargs['y_col'], list):
        y = None
    else:
        y = data_frame[kwargs['y_col']].ravel()

    # load groups
    if groups:
        groups_selector = (params['experiment_schemes']['test_split']
                                 ['split_algos']).pop('groups_selector')

        header = 'infer' if groups_selector['header_g'] else None
        column_option = \
            (groups_selector['column_selector_options_g']
                            ['selected_column_selector_option_g'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = groups_selector['column_selector_options_g']['col_g']
        else:
            c = None

        if groups == infile_dataframe:
            groups = data_frame

        groups = read_columns(
                groups,
                c=c,
                c_option=column_option,
                sep='\t',
                header=header,
                parse_dates=True)
        groups = groups.ravel()

    # handle test (first) split
    test_split_options = (params['experiment_schemes']
                                ['test_split']['split_algos'])

    if test_split_options['shuffle'] == 'group':
        test_split_options['labels'] = groups
    if test_split_options['shuffle'] == 'stratified':
        if y is not None:
            test_split_options['labels'] = y
        else:
            raise ValueError("Stratified shuffle split is not "
                             "applicable on empty target values or "
                             "multiple output targets!")

    X_train, X_test, y_train, y_test, groups_train, groups_test = \
        train_test_split_none(X, y, groups, **test_split_options)

    exp_scheme = params['experiment_schemes']['selected_exp_scheme']

    # handle validation (second) split
    if exp_scheme == 'train_val_test':
        val_split_options = (params['experiment_schemes']
                                   ['val_split']['split_algos'])

        if val_split_options['shuffle'] == 'group':
            val_split_options['labels'] = groups_train
        if val_split_options['shuffle'] == 'stratified':
            if y_train is not None:
                val_split_options['labels'] = y_train
            else:
                raise ValueError("Stratified shuffle split is not "
                                 "applicable on empty target values!")

        X_train, X_val, y_train, y_val, groups_train, groups_val = \
            train_test_split_none(X_train, y_train, groups_train,
                                  **val_split_options)

        # In image data generator, `y_val` must be None
        # labels will be retrived in generator.
        estimator.fit(X_train, y_train,
                      validation_data=(X_val, ))

    else:
        estimator.fit(X_train, y_train,
                      validation_data=(X_test, ))

    # evaluation
    scores = {}
    steps = estimator.prediction_steps
    batch_size = estimator.batch_size
    generator = estimator.data_generator_.flow(X_test, y=y_test,
                                               batch_size=batch_size)

    # keras metrics evaluation
    # handle scorer, convert to scorer dict
    generator.reset()
    score_results = estimator.model_.evaluate_generator(generator,
                                                        steps=steps)
    metrics_names = estimator.model_.metrics_names
    if not isinstance(metrics_names, list):
        scores[metrics_names] = score_results
    else:
        for i, mm in enumerate(metrics_names):
            scores[mm] = score_results[i]

    # for sklearn metrics
    scoring = params['experiment_schemes']['metrics']['scoring']
    if scoring['primary_scoring'] != 'default' or outfile_y_true:
        generator.reset()
        predictions, y_true = _predict_generator(estimator.model_,
                                                 generator,
                                                 steps=steps)

    if scoring['primary_scoring'] != 'default':
        scorer = get_scoring(scoring)
        scorer, _ = _check_multimetric_scoring(estimator,
                                               scoring=scorer)
        sk_scores = gen_compute_scores(y_true, predictions, scorer,
                                       is_multimetric=True)
        scores.update(sk_scores)

    if outfile_y_true:
        try:
            pd.DataFrame(y_true).to_csv(outfile_y_true, sep='\t',
                                        index=False)
            pd.DataFrame(predictions).astype(np.float32).to_csv(
                outfile_y_preds, sep='\t', index=False,
                float_format='%g', chunksize=10000)
        except Exception as e:
            print("Error in saving predictions: %s" % e)

    # handle output
    for name, score in scores.items():
        scores[name] = [score]
    df = pd.DataFrame(scores)
    df = df[sorted(df.columns)]
    df.to_csv(path_or_buf=outfile_result, sep='\t', header=True,
              index=False)

    if outfile_object:
        main_est = get_main_estimator(estimator)

        if hasattr(main_est, 'model_') \
                and hasattr(main_est, 'save_weights'):
            if outfile_weights:
                main_est.save_weights(outfile_weights)
            del main_est.model_
            del main_est.fit_params
            del main_est.model_class_
            if getattr(main_est, 'data_generator_', None):
                del main_est.data_generator_

        with open(outfile_object, 'wb') as output_handler:
            pickle.dump(estimator, output_handler,
                        pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--estimator", dest="infile_estimator")
    aparser.add_argument("-X", "--infile_images", dest="infile_images")
    aparser.add_argument("-y", "--infile_dataframe", dest="infile_dataframe")
    aparser.add_argument("-O", "--outfile_result", dest="outfile_result")
    aparser.add_argument("-o", "--outfile_object", dest="outfile_object")
    aparser.add_argument("-w", "--outfile_weights", dest="outfile_weights")
    aparser.add_argument("-l", "--outfile_y_true", dest="outfile_y_true")
    aparser.add_argument("-p", "--outfile_y_preds", dest="outfile_y_preds")
    aparser.add_argument("-g", "--groups", dest="groups")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile_images,
         args.infile_dataframe, args.outfile_result,
         outfile_object=args.outfile_object,
         outfile_weights=args.outfile_weights,
         outfile_y_true=args.outfile_y_true,
         outfile_y_preds=args.outfile_y_preds,
         groups=args.groups)
