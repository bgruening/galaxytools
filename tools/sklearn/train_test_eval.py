import argparse
import joblib
import json
import numpy as np
import pandas as pd
import pickle
import warnings
from itertools import chain
from scipy.io import mmread
from sklearn.base import clone
from sklearn import (cluster, compose, decomposition, ensemble,
                     feature_extraction, feature_selection,
                     gaussian_process, kernel_approximation, metrics,
                     model_selection, naive_bayes, neighbors,
                     pipeline, preprocessing, svm, linear_model,
                     tree, discriminant_analysis)
from sklearn.exceptions import FitFailedWarning
from sklearn.metrics.scorer import _check_multimetric_scoring
from sklearn.model_selection._validation import _score, cross_validate
from sklearn.model_selection import _search, _validation
from sklearn.utils import indexable, safe_indexing

from galaxy_ml.model_validations import train_test_split
from galaxy_ml.utils import (SafeEval, get_scoring, load_model,
                             read_columns, try_get_attr, get_module)


_fit_and_score = try_get_attr('galaxy_ml.model_validations', '_fit_and_score')
setattr(_search, '_fit_and_score', _fit_and_score)
setattr(_validation, '_fit_and_score', _fit_and_score)

N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))
CACHE_DIR = './cached'
NON_SEARCHABLE = ('n_jobs', 'pre_dispatch', 'memory', '_path',
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


def main(inputs, infile_estimator, infile1, infile2,
         outfile_result, outfile_object=None,
         outfile_weights=None, groups=None,
         ref_seq=None, intervals=None, targets=None,
         fasta_path=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_estimator : str
        File path to estimator

    infile1 : str
        File path to dataset containing features

    infile2 : str
        File path to dataset containing target values

    outfile_result : str
        File path to save the results, either cv_results or test result

    outfile_object : str, optional
        File path to save searchCV object

    outfile_weights : str, optional
        File path to save deep learning model weights

    groups : str
        File path to dataset containing groups labels

    ref_seq : str
        File path to dataset containing genome sequence file

    intervals : str
        File path to dataset containing interval file

    targets : str
        File path to dataset compressed target bed file

    fasta_path : str
        File path to dataset containing fasta file
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    #  load estimator
    with open(infile_estimator, 'rb') as estimator_handler:
        estimator = load_model(estimator_handler)

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
            parse_dates=True)
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
                parse_dates=True)
        groups = groups.ravel()

    # del loaded_df
    del loaded_df

    # handle memory
    memory = joblib.Memory(location=CACHE_DIR, verbose=0)
    # cache iraps_core fits could increase search speed significantly
    if estimator.__class__.__name__ == 'IRAPSClassifier':
        estimator.set_params(memory=memory)
    else:
        # For iraps buried in pipeline
        new_params = {}
        for p, v in estimator_params.items():
            if p.endswith('memory'):
                # for case of `__irapsclassifier__memory`
                if len(p) > 8 and p[:-8].endswith('irapsclassifier'):
                    # cache iraps_core fits could increase search
                    # speed significantly
                    new_params[p] = memory
                # security reason, we don't want memory being
                # modified unexpectedly
                elif v:
                    new_params[p] = None
            # handle n_jobs
            elif p.endswith('n_jobs'):
                # For now, 1 CPU is suggested for iprasclassifier
                if len(p) > 8 and p[:-8].endswith('irapsclassifier'):
                    new_params[p] = 1
                else:
                    new_params[p] = N_JOBS
            # for security reason, types of callback are limited
            elif p.endswith('callbacks'):
                for cb in v:
                    cb_type = cb['callback_selection']['callback_type']
                    if cb_type not in ALLOWED_CALLBACKS:
                        raise ValueError(
                            "Prohibited callback type: %s!" % cb_type)

        estimator.set_params(**new_params)

    # handle scorer, convert to scorer dict
    scoring = params['experiment_schemes']['metrics']['scoring']
    scorer = get_scoring(scoring)
    scorer, _ = _check_multimetric_scoring(estimator, scoring=scorer)

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
    if hasattr(estimator, 'validation_data'):
        if exp_scheme == 'train_val_test':
            estimator.fit(X_train, y_train,
                          validation_data=(X_val, y_val))
        else:
            estimator.fit(X_train, y_train,
                          validation_data=(X_test, y_test))
    else:
        estimator.fit(X_train, y_train)

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
    df.to_csv(path_or_buf=outfile_result, sep='\t',
              header=True, index=False)

    memory.clear(warn=False)

    if outfile_object:
        main_est = estimator
        if isinstance(estimator, pipeline.Pipeline):
            main_est = estimator.steps[-1][-1]

        if hasattr(main_est, 'model_') \
                and hasattr(main_est, 'save_weights'):
            if outfile_weights:
                main_est.save_weights(outfile_weights)
            del main_est.model_

        with open(outfile_object, 'wb') as output_handler:
            pickle.dump(estimator, output_handler,
                        pickle.HIGHEST_PROTOCOL)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-e", "--estimator", dest="infile_estimator")
    aparser.add_argument("-X", "--infile1", dest="infile1")
    aparser.add_argument("-y", "--infile2", dest="infile2")
    aparser.add_argument("-O", "--outfile_result", dest="outfile_result")
    aparser.add_argument("-o", "--outfile_object", dest="outfile_object")
    aparser.add_argument("-w", "--outfile_weights", dest="outfile_weights")
    aparser.add_argument("-g", "--groups", dest="groups")
    aparser.add_argument("-r", "--ref_seq", dest="ref_seq")
    aparser.add_argument("-b", "--intervals", dest="intervals")
    aparser.add_argument("-t", "--targets", dest="targets")
    aparser.add_argument("-f", "--fasta_path", dest="fasta_path")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile1, args.infile2,
         args.outfile_result, outfile_object=args.outfile_object,
         outfile_weights=args.outfile_weights, groups=args.groups,
         ref_seq=args.ref_seq, intervals=args.intervals,
         targets=args.targets, fasta_path=args.fasta_path)
