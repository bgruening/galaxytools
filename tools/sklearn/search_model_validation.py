import argparse
import collections
import imblearn
import joblib
import json
import numpy as np
import os
import pandas as pd
import pickle
import skrebate
import sys
import warnings
from scipy.io import mmread
from sklearn import (cluster, decomposition, feature_selection,
                     kernel_approximation, model_selection, preprocessing)
from sklearn.exceptions import FitFailedWarning
from sklearn.model_selection._validation import _score, cross_validate
from sklearn.model_selection import _search, _validation
from sklearn.pipeline import Pipeline

from galaxy_ml.binarize_target import IRAPSClassifier
from galaxy_ml.utils import (SafeEval, get_cv, get_scoring, load_model,
                             read_columns, try_get_attr, get_module,
                             clean_params, get_main_estimator)


N_JOBS = int(os.environ.get('GALAXY_SLOTS', 1))
# handle  disk cache
CACHE_DIR = os.path.join(os.getcwd(), 'cached')
del os
NON_SEARCHABLE = ('n_jobs', 'pre_dispatch', 'memory', '_path', '_dir',
                  'nthread', 'callbacks')


def _eval_search_params(params_builder):
    search_params = {}

    for p in params_builder['param_set']:
        search_list = p['sp_list'].strip()
        if search_list == '':
            continue

        param_name = p['sp_name']
        if param_name.lower().endswith(NON_SEARCHABLE):
            print("Warning: `%s` is not eligible for search and was "
                  "omitted!" % param_name)
            continue

        if not search_list.startswith(':'):
            safe_eval = SafeEval(load_scipy=True, load_numpy=True)
            ev = safe_eval(search_list)
            search_params[param_name] = ev
        else:
            # Have `:` before search list, asks for estimator evaluatio
            safe_eval_es = SafeEval(load_estimators=True)
            search_list = search_list[1:].strip()
            # TODO maybe add regular express check
            ev = safe_eval_es(search_list)
            preprocessings = (
                preprocessing.StandardScaler(), preprocessing.Binarizer(),
                preprocessing.MaxAbsScaler(),
                preprocessing.Normalizer(), preprocessing.MinMaxScaler(),
                preprocessing.PolynomialFeatures(),
                preprocessing.RobustScaler(), feature_selection.SelectKBest(),
                feature_selection.GenericUnivariateSelect(),
                feature_selection.SelectPercentile(),
                feature_selection.SelectFpr(), feature_selection.SelectFdr(),
                feature_selection.SelectFwe(),
                feature_selection.VarianceThreshold(),
                decomposition.FactorAnalysis(random_state=0),
                decomposition.FastICA(random_state=0),
                decomposition.IncrementalPCA(),
                decomposition.KernelPCA(random_state=0, n_jobs=N_JOBS),
                decomposition.LatentDirichletAllocation(
                    random_state=0, n_jobs=N_JOBS),
                decomposition.MiniBatchDictionaryLearning(
                    random_state=0, n_jobs=N_JOBS),
                decomposition.MiniBatchSparsePCA(
                    random_state=0, n_jobs=N_JOBS),
                decomposition.NMF(random_state=0),
                decomposition.PCA(random_state=0),
                decomposition.SparsePCA(random_state=0, n_jobs=N_JOBS),
                decomposition.TruncatedSVD(random_state=0),
                kernel_approximation.Nystroem(random_state=0),
                kernel_approximation.RBFSampler(random_state=0),
                kernel_approximation.AdditiveChi2Sampler(),
                kernel_approximation.SkewedChi2Sampler(random_state=0),
                cluster.FeatureAgglomeration(),
                skrebate.ReliefF(n_jobs=N_JOBS),
                skrebate.SURF(n_jobs=N_JOBS),
                skrebate.SURFstar(n_jobs=N_JOBS),
                skrebate.MultiSURF(n_jobs=N_JOBS),
                skrebate.MultiSURFstar(n_jobs=N_JOBS),
                imblearn.under_sampling.ClusterCentroids(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.CondensedNearestNeighbour(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.EditedNearestNeighbours(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.RepeatedEditedNearestNeighbours(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.AllKNN(random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.InstanceHardnessThreshold(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.NearMiss(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.NeighbourhoodCleaningRule(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.OneSidedSelection(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.under_sampling.RandomUnderSampler(
                    random_state=0),
                imblearn.under_sampling.TomekLinks(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.over_sampling.ADASYN(random_state=0, n_jobs=N_JOBS),
                imblearn.over_sampling.RandomOverSampler(random_state=0),
                imblearn.over_sampling.SMOTE(random_state=0, n_jobs=N_JOBS),
                imblearn.over_sampling.SVMSMOTE(random_state=0, n_jobs=N_JOBS),
                imblearn.over_sampling.BorderlineSMOTE(
                    random_state=0, n_jobs=N_JOBS),
                imblearn.over_sampling.SMOTENC(
                    categorical_features=[], random_state=0, n_jobs=N_JOBS),
                imblearn.combine.SMOTEENN(random_state=0),
                imblearn.combine.SMOTETomek(random_state=0))
            newlist = []
            for obj in ev:
                if obj is None:
                    newlist.append(None)
                elif obj == 'all_0':
                    newlist.extend(preprocessings[0:35])
                elif obj == 'sk_prep_all':      # no KernalCenter()
                    newlist.extend(preprocessings[0:7])
                elif obj == 'fs_all':
                    newlist.extend(preprocessings[7:14])
                elif obj == 'decomp_all':
                    newlist.extend(preprocessings[14:25])
                elif obj == 'k_appr_all':
                    newlist.extend(preprocessings[25:29])
                elif obj == 'reb_all':
                    newlist.extend(preprocessings[30:35])
                elif obj == 'imb_all':
                    newlist.extend(preprocessings[35:54])
                elif type(obj) is int and -1 < obj < len(preprocessings):
                    newlist.append(preprocessings[obj])
                elif hasattr(obj, 'get_params'):       # user uploaded object
                    if 'n_jobs' in obj.get_params():
                        newlist.append(obj.set_params(n_jobs=N_JOBS))
                    else:
                        newlist.append(obj)
                else:
                    sys.exit("Unsupported estimator type: %r" % (obj))

            search_params[param_name] = newlist

    return search_params


def _handle_X_y(estimator, params, infile1, infile2, loaded_df={},
                ref_seq=None, intervals=None, targets=None,
                fasta_path=None):
    """read inputs

    Params
    -------
    estimator : estimator object
    params : dict
        Galaxy tool parameter inputs
    infile1 : str
        File path to dataset containing features
    infile2 : str
        File path to dataset containing target values
    loaded_df : dict
        Contains loaded DataFrame objects with file path as keys
    ref_seq : str
        File path to dataset containing genome sequence file
    interval : str
        File path to dataset containing interval file
    targets : str
        File path to dataset compressed target bed file
    fasta_path : str
        File path to dataset containing fasta file


    Returns
    -------
    estimator : estimator object after setting new attributes
    X : numpy array
    y : numpy array
    """
    estimator_params = estimator.get_params()

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

        if df_key in loaded_df:
            infile1 = loaded_df[df_key]

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

    return estimator, X, y


def _do_outer_cv(searcher, X, y, outer_cv, scoring, error_score='raise',
                 outfile=None):
    """Do outer cross-validation for nested CV

    Parameters
    ----------
    searcher : object
        SearchCV object
    X : numpy array
        Containing features
    y : numpy array
        Target values or labels
    outer_cv : int or CV splitter
        Control the cv splitting
    scoring : object
        Scorer
    error_score: str, float or numpy float
        Whether to raise fit error or return an value
    outfile : str
        File path to store the restuls
    """
    if error_score == 'raise':
        rval = cross_validate(
            searcher, X, y, scoring=scoring,
            cv=outer_cv, n_jobs=N_JOBS, verbose=0,
            error_score=error_score)
    else:
        warnings.simplefilter('always', FitFailedWarning)
        with warnings.catch_warnings(record=True) as w:
            try:
                rval = cross_validate(
                    searcher, X, y,
                    scoring=scoring,
                    cv=outer_cv, n_jobs=N_JOBS,
                    verbose=0,
                    error_score=error_score)
            except ValueError:
                pass
            for warning in w:
                print(repr(warning.message))

    keys = list(rval.keys())
    for k in keys:
        if k.startswith('test'):
            rval['mean_' + k] = np.mean(rval[k])
            rval['std_' + k] = np.std(rval[k])
        if k.endswith('time'):
            rval.pop(k)
    rval = pd.DataFrame(rval)
    rval = rval[sorted(rval.columns)]
    rval.to_csv(path_or_buf=outfile, sep='\t', header=True, index=False)


def _do_train_test_split_val(searcher, X, y, params, error_score='raise',
                             primary_scoring=None, groups=None,
                             outfile=None):
    """ do train test split, searchCV validates on the train and then use
    the best_estimator_ to evaluate on the test

    Returns
    --------
    Fitted SearchCV object
    """
    train_test_split = try_get_attr(
        'galaxy_ml.model_validations', 'train_test_split')
    split_options = params['outer_split']

    # splits
    if split_options['shuffle'] == 'stratified':
        split_options['labels'] = y
        X, X_test, y, y_test = train_test_split(X, y, **split_options)
    elif split_options['shuffle'] == 'group':
        if groups is None:
            raise ValueError("No group based CV option was choosen for "
                             "group shuffle!")
        split_options['labels'] = groups
        if y is None:
            X, X_test, groups, _ =\
                train_test_split(X, groups, **split_options)
        else:
            X, X_test, y, y_test, groups, _ =\
                train_test_split(X, y, groups, **split_options)
    else:
        if split_options['shuffle'] == 'None':
            split_options['shuffle'] = None
        X, X_test, y, y_test =\
            train_test_split(X, y, **split_options)

    if error_score == 'raise':
        searcher.fit(X, y, groups=groups)
    else:
        warnings.simplefilter('always', FitFailedWarning)
        with warnings.catch_warnings(record=True) as w:
            try:
                searcher.fit(X, y, groups=groups)
            except ValueError:
                pass
            for warning in w:
                print(repr(warning.message))

    scorer_ = searcher.scorer_
    if isinstance(scorer_, collections.Mapping):
        is_multimetric = True
    else:
        is_multimetric = False

    best_estimator_ = getattr(searcher, 'best_estimator_')

    # TODO Solve deep learning models in pipeline
    if best_estimator_.__class__.__name__ == 'KerasGBatchClassifier':
        test_score = best_estimator_.evaluate(
            X_test, scorer=scorer_, is_multimetric=is_multimetric)
    else:
        test_score = _score(best_estimator_, X_test,
                            y_test, scorer_,
                            is_multimetric=is_multimetric)

    if not is_multimetric:
        test_score = {primary_scoring: test_score}
    for key, value in test_score.items():
        test_score[key] = [value]
    result_df = pd.DataFrame(test_score)
    result_df.to_csv(path_or_buf=outfile, sep='\t', header=True,
                     index=False)

    return searcher


def _set_memory(estimator, memory):
    """set memeory cache

    Parameters
    ----------
    estimator : python object
    memory : joblib.Memory object

    Returns
    -------
    estimator : estimator object after setting new attributes
    """
    if isinstance(estimator, IRAPSClassifier):
        estimator.set_params(memory=memory)
        return estimator

    estimator_params = estimator.get_params()

    new_params = {}
    for k in estimator_params.keys():
        if k.endswith('irapsclassifier__memory'):
            new_params[k] = memory

    estimator.set_params(**new_params)

    return estimator


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
        File path to save model weights

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

    # store read dataframe object
    loaded_df = {}

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    # Override the refit parameter
    params['search_schemes']['options']['refit'] = True \
        if (params['save'] != 'nope' or
            params['outer_split']['split_mode'] == 'nested_cv') else False

    with open(infile_estimator, 'rb') as estimator_handler:
        estimator = load_model(estimator_handler)

    if estimator.__class__.__name__ == 'KerasGBatchClassifier':
        _fit_and_score = try_get_attr('galaxy_ml.model_validations',
                                      '_fit_and_score')

        setattr(_search, '_fit_and_score', _fit_and_score)
        setattr(_validation, '_fit_and_score', _fit_and_score)

    optimizer = params['search_schemes']['selected_search_scheme']
    optimizer = getattr(model_selection, optimizer)

    # handle gridsearchcv options
    options = params['search_schemes']['options']

    if groups:
        header = 'infer' if (options['cv_selector']['groups_selector']
                                    ['header_g']) else None
        column_option = (options['cv_selector']['groups_selector']
                                ['column_selector_options_g']
                                ['selected_column_selector_option_g'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = (options['cv_selector']['groups_selector']
                        ['column_selector_options_g']['col_g'])
        else:
            c = None

        df_key = groups + repr(header)

        groups = pd.read_csv(groups, sep='\t', header=header,
                             parse_dates=True)
        loaded_df[df_key] = groups

        groups = read_columns(
                groups,
                c=c,
                c_option=column_option,
                sep='\t',
                header=header,
                parse_dates=True)
        groups = groups.ravel()
        options['cv_selector']['groups_selector'] = groups

    splitter, groups = get_cv(options.pop('cv_selector'))
    options['cv'] = splitter
    primary_scoring = options['scoring']['primary_scoring']
    options['scoring'] = get_scoring(options['scoring'])
    if options['error_score']:
        options['error_score'] = 'raise'
    else:
        options['error_score'] = np.NaN
    if options['refit'] and isinstance(options['scoring'], dict):
        options['refit'] = primary_scoring
    if 'pre_dispatch' in options and options['pre_dispatch'] == '':
        options['pre_dispatch'] = None

    params_builder = params['search_schemes']['search_params_builder']
    param_grid = _eval_search_params(params_builder)

    estimator = clean_params(estimator)

    # save the SearchCV object without fit
    if params['save'] == 'save_no_fit':
        searcher = optimizer(estimator, param_grid, **options)
        print(searcher)
        with open(outfile_object, 'wb') as output_handler:
            pickle.dump(searcher, output_handler,
                        pickle.HIGHEST_PROTOCOL)
        return 0

    # read inputs and loads new attributes, like paths
    estimator, X, y = _handle_X_y(estimator, params, infile1, infile2,
                                  loaded_df=loaded_df, ref_seq=ref_seq,
                                  intervals=intervals, targets=targets,
                                  fasta_path=fasta_path)

    # cache iraps_core fits could increase search speed significantly
    memory = joblib.Memory(location=CACHE_DIR, verbose=0)
    estimator = _set_memory(estimator, memory)

    searcher = optimizer(estimator, param_grid, **options)

    split_mode = params['outer_split'].pop('split_mode')

    if split_mode == 'nested_cv':
        outer_cv, _ = get_cv(params['outer_split']['cv_selector'])
        # nested CV, outer cv using cross_validate
        if options['error_score'] == 'raise':
            rval = cross_validate(
                searcher, X, y, scoring=options['scoring'],
                cv=outer_cv, n_jobs=N_JOBS,
                verbose=options['verbose'],
                return_estimator=(params['save'] == 'save_estimator'),
                error_score=options['error_score'],
                return_train_score=True)
        else:
            warnings.simplefilter('always', FitFailedWarning)
            with warnings.catch_warnings(record=True) as w:
                try:
                    rval = cross_validate(
                        searcher, X, y,
                        scoring=options['scoring'],
                        cv=outer_cv, n_jobs=N_JOBS,
                        verbose=options['verbose'],
                        return_estimator=(params['save'] == 'save_estimator'),
                        error_score=options['error_score'],
                        return_train_score=True)
                except ValueError:
                    pass
                for warning in w:
                    print(repr(warning.message))

        fitted_searchers = rval.pop('estimator', [])
        if fitted_searchers:
            import os
            pwd = os.getcwd()
            save_dir = os.path.join(pwd, 'cv_results_in_folds')
            try:
                os.mkdir(save_dir)
                for idx, obj in enumerate(fitted_searchers):
                    target_name = 'cv_results_' + '_' + 'split%d' % idx
                    target_path = os.path.join(pwd, save_dir, target_name)
                    cv_results_ = getattr(obj, 'cv_results_', None)
                    if not cv_results_:
                        print("%s is not available" % target_name)
                        continue
                    cv_results_ = pd.DataFrame(cv_results_)
                    cv_results_ = cv_results_[sorted(cv_results_.columns)]
                    cv_results_.to_csv(target_path, sep='\t', header=True,
                                       index=False)
            except Exception as e:
                print(e)
            finally:
                del os

        keys = list(rval.keys())
        for k in keys:
            if k.startswith('test'):
                rval['mean_' + k] = np.mean(rval[k])
                rval['std_' + k] = np.std(rval[k])
            if k.endswith('time'):
                rval.pop(k)
        rval = pd.DataFrame(rval)
        rval = rval[sorted(rval.columns)]
        rval.to_csv(path_or_buf=outfile_result, sep='\t', header=True,
                    index=False)

        return 0

        # deprecate train test split mode
        """searcher = _do_train_test_split_val(
            searcher, X, y, params,
            primary_scoring=primary_scoring,
            error_score=options['error_score'],
            groups=groups,
            outfile=outfile_result)"""

    # no outer split
    else:
        searcher.set_params(n_jobs=N_JOBS)
        if options['error_score'] == 'raise':
            searcher.fit(X, y, groups=groups)
        else:
            warnings.simplefilter('always', FitFailedWarning)
            with warnings.catch_warnings(record=True) as w:
                try:
                    searcher.fit(X, y, groups=groups)
                except ValueError:
                    pass
                for warning in w:
                    print(repr(warning.message))

        cv_results = pd.DataFrame(searcher.cv_results_)
        cv_results = cv_results[sorted(cv_results.columns)]
        cv_results.to_csv(path_or_buf=outfile_result, sep='\t',
                          header=True, index=False)

    memory.clear(warn=False)

    # output best estimator, and weights if applicable
    if outfile_object:
        best_estimator_ = getattr(searcher, 'best_estimator_', None)
        if not best_estimator_:
            warnings.warn("GridSearchCV object has no attribute "
                          "'best_estimator_', because either it's "
                          "nested gridsearch or `refit` is False!")
            return

        # clean prams
        best_estimator_ = clean_params(best_estimator_)

        main_est = get_main_estimator(best_estimator_)

        if hasattr(main_est, 'model_') \
                and hasattr(main_est, 'save_weights'):
            if outfile_weights:
                main_est.save_weights(outfile_weights)
            del main_est.model_
            del main_est.fit_params
            del main_est.model_class_
            main_est.callbacks = []
            if getattr(main_est, 'data_generator_', None):
                del main_est.data_generator_

        with open(outfile_object, 'wb') as output_handler:
            print("Best estimator is saved: %s " % repr(best_estimator_))
            pickle.dump(best_estimator_, output_handler,
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
