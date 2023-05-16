import argparse
import json
import os
import warnings
from itertools import chain

from galaxy_ml.keras_galaxy_models import (
    KerasGBatchClassifier, _predict_generator,
)
from galaxy_ml.model_persist import dump_model_to_h5, load_model_from_h5
from galaxy_ml.model_validations import train_test_split
from galaxy_ml.utils import (
    SafeEval, clean_params, gen_compute_scores, get_main_estimator,
    get_module, get_scoring, read_columns,
)

import joblib

import numpy as np

import pandas as pd

from scipy.io import mmread

from sklearn.metrics._scorer import _check_multimetric_scoring
from sklearn.model_selection._validation import _score
from sklearn.pipeline import Pipeline
from sklearn.utils import _safe_indexing, indexable


N_JOBS = int(os.environ.get('GALAXY_SLOTS', 1))
CACHE_DIR = os.path.join(os.getcwd(), 'cached')
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
            (_safe_indexing(a, train),
             _safe_indexing(a, test)) for a in new_arrays))
    else:
        rval = train_test_split(*new_arrays, **kwargs)

    for pos in nones:
        rval[pos * 2: 2] = [None, None]

    return rval


def _evaluate_keras_and_sklearn_scores(estimator, data_generator, X,
                                       y=None, sk_scoring=None,
                                       steps=None, batch_size=32,
                                       return_predictions=False):
    """output scores for bother keras and sklearn metrics

    Parameters
    -----------
    estimator : object
        Fitted `galaxy_ml.keras_galaxy_models.KerasGBatchClassifier`.
    data_generator : object
        From `galaxy_ml.preprocessors.ImageDataFrameBatchGenerator`.
    X : 2-D array
        Contains indecies of images that need to be evaluated.
    y : None
        Target value.
    sk_scoring : dict
        Galaxy tool input parameters.
    steps : integer or None
        Evaluation/prediction steps before stop.
    batch_size : integer
        Number of samples in a batch
    return_predictions : bool, default is False
        Whether to return predictions and true labels.
    """
    scores = {}

    generator = data_generator.flow(X, y=y, batch_size=batch_size)
    # keras metrics evaluation
    # handle scorer, convert to scorer dict
    generator.reset()
    score_results = estimator.model_.evaluate_generator(generator,
                                                        steps=steps)
    metrics_names = estimator.model_.metrics_names
    if not isinstance(metrics_names, list):
        scores[metrics_names] = score_results
    else:
        scores = dict(zip(metrics_names, score_results))

    if sk_scoring['primary_scoring'] == 'default' and\
            not return_predictions:
        return scores

    generator.reset()
    predictions, y_true = _predict_generator(estimator.model_,
                                             generator,
                                             steps=steps)

    # for sklearn metrics
    if sk_scoring['primary_scoring'] != 'default':
        scorer = get_scoring(sk_scoring)
        if not isinstance(scorer, (dict, list)):
            scorer = [sk_scoring['primary_scoring']]
        scorer = _check_multimetric_scoring(estimator, scoring=scorer)
        sk_scores = gen_compute_scores(y_true, predictions, scorer)
        scores.update(sk_scores)

    if return_predictions:
        return scores, predictions, y_true
    else:
        return scores, None, None


def main(inputs, infile_estimator, infile1, infile2,
         outfile_result, outfile_object=None,
         outfile_y_true=None,
         outfile_y_preds=None, groups=None,
         ref_seq=None, intervals=None, targets=None,
         fasta_path=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter.

    infile_estimator : str
        File path to estimator.

    infile1 : str
        File path to dataset containing features.

    infile2 : str
        File path to dataset containing target values.

    outfile_result : str
        File path to save the results, either cv_results or test result.

    outfile_object : str, optional
        File path to save searchCV object.

    outfile_y_true : str, optional
        File path to target values for prediction.

    outfile_y_preds : str, optional
        File path to save predictions.

    groups : str
        File path to dataset containing groups labels.

    ref_seq : str
        File path to dataset containing genome sequence file.

    intervals : str
        File path to dataset containing interval file.

    targets : str
        File path to dataset compressed target bed file.

    fasta_path : str
        File path to dataset containing fasta file.
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    #  load estimator
    estimator = load_model_from_h5(infile_estimator)

    estimator = clean_params(estimator)

    # swap hyperparameter
    swapping = params['experiment_schemes']['hyperparams_swapping']
    swap_params = _eval_swap_params(swapping)
    estimator.set_params(**swap_params)

    estimator_params = estimator.get_params()

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

    # fasta_file input
    elif input_type == 'seq_fasta':
        pyfaidx = get_module('pyfaidx')
        sequences = pyfaidx.Fasta(fasta_path)
        n_seqs = len(sequences.keys())
        X = np.arange(n_seqs)[:, np.newaxis]
        for param in estimator_params.keys():
            if param.endswith('fasta_path'):
                estimator.set_params(
                    **{param: fasta_path})
                break
        else:
            raise ValueError(
                "The selected estimator doesn't support "
                "fasta file input! Please consider using "
                "KerasGBatchClassifier with "
                "FastaDNABatchGenerator/FastaProteinBatchGenerator "
                "or having GenomeOneHotEncoder/ProteinOneHotEncoder "
                "in pipeline!")

    elif input_type == 'refseq_and_interval':
        path_params = {
            'data_batch_generator__ref_genome_path': ref_seq,
            'data_batch_generator__intervals_path': intervals,
            'data_batch_generator__target_path': targets
        }
        estimator.set_params(**path_params)
        n_intervals = sum(1 for line in open(intervals))
        X = np.arange(n_intervals)[:, np.newaxis]

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
        parse_dates=True,
    )
    if len(y.shape) == 2 and y.shape[1] == 1:
        y = y.ravel()
    if input_type == 'refseq_and_interval':
        estimator.set_params(
            data_batch_generator__features=y.ravel().tolist())
        y = None
    # end y

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

        df_key = groups + repr(header)
        if df_key in loaded_df:
            groups = loaded_df[df_key]

        groups = read_columns(
            groups,
            c=c,
            c_option=column_option,
            sep='\t',
            header=header,
            parse_dates=True,
        )
        groups = groups.ravel()

    # del loaded_df
    del loaded_df

    # cache iraps_core fits could increase search speed significantly
    memory = joblib.Memory(location=CACHE_DIR, verbose=0)
    main_est = get_main_estimator(estimator)
    if main_est.__class__.__name__ == 'IRAPSClassifier':
        main_est.set_params(memory=memory)

    # handle scorer, convert to scorer dict
    scoring = params['experiment_schemes']['metrics']['scoring']
    scorer = get_scoring(scoring)
    if not isinstance(scorer, (dict, list)):
        scorer = [scoring['primary_scoring']]
    scorer = _check_multimetric_scoring(estimator, scoring=scorer)

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
                             "applicable on empty target values!")

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

    # train and eval
    if hasattr(estimator, 'config') and hasattr(estimator, 'model_type'):
        if exp_scheme == 'train_val_test':
            estimator.fit(X_train, y_train,
                          validation_data=(X_val, y_val))
        else:
            estimator.fit(X_train, y_train,
                          validation_data=(X_test, y_test))
    else:
        estimator.fit(X_train, y_train)

    if isinstance(estimator, KerasGBatchClassifier):
        scores = {}
        steps = estimator.prediction_steps
        batch_size = estimator.batch_size
        data_generator = estimator.data_generator_

        scores, predictions, y_true = _evaluate_keras_and_sklearn_scores(
            estimator, data_generator, X_test, y=y_test,
            sk_scoring=scoring, steps=steps, batch_size=batch_size,
            return_predictions=bool(outfile_y_true))

    else:
        scores = {}
        if hasattr(estimator, 'model_') \
                and hasattr(estimator.model_, 'metrics_names'):
            batch_size = estimator.batch_size
            score_results = estimator.model_.evaluate(X_test, y=y_test,
                                                      batch_size=batch_size,
                                                      verbose=0)
            metrics_names = estimator.model_.metrics_names
            if not isinstance(metrics_names, list):
                scores[metrics_names] = score_results
            else:
                scores = dict(zip(metrics_names, score_results))

        if hasattr(estimator, 'predict_proba'):
            predictions = estimator.predict_proba(X_test)
        else:
            predictions = estimator.predict(X_test)

        y_true = y_test
        sk_scores = _score(estimator, X_test, y_test, scorer)
        scores.update(sk_scores)

    # handle output
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
    df.to_csv(path_or_buf=outfile_result, sep='\t',
              header=True, index=False)

    memory.clear(warn=False)

    if outfile_object:
        dump_model_to_h5(estimator, outfile_object)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--estimator", dest="infile_estimator")
    aparser.add_argument("-X", "--infile1", dest="infile1")
    aparser.add_argument("-y", "--infile2", dest="infile2")
    aparser.add_argument("-O", "--outfile_result", dest="outfile_result")
    aparser.add_argument("-o", "--outfile_object", dest="outfile_object")
    aparser.add_argument("-l", "--outfile_y_true", dest="outfile_y_true")
    aparser.add_argument("-p", "--outfile_y_preds", dest="outfile_y_preds")
    aparser.add_argument("-g", "--groups", dest="groups")
    aparser.add_argument("-r", "--ref_seq", dest="ref_seq")
    aparser.add_argument("-b", "--intervals", dest="intervals")
    aparser.add_argument("-t", "--targets", dest="targets")
    aparser.add_argument("-f", "--fasta_path", dest="fasta_path")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile1, args.infile2,
         args.outfile_result, outfile_object=args.outfile_object,
         outfile_y_true=args.outfile_y_true,
         outfile_y_preds=args.outfile_y_preds,
         groups=args.groups,
         ref_seq=args.ref_seq, intervals=args.intervals,
         targets=args.targets, fasta_path=args.fasta_path)
