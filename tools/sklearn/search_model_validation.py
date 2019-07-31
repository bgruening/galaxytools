import argparse
import collections
import imblearn
import joblib
import json
import numpy as np
import pandas as pd
import pickle
import skrebate
import sklearn
import sys
import xgboost
import warnings
from imblearn import under_sampling, over_sampling, combine
from scipy.io import mmread
from mlxtend import classifier, regressor
from sklearn.base import clone
from sklearn import (cluster, compose, decomposition, ensemble,
                     feature_extraction, feature_selection,
                     gaussian_process, kernel_approximation, metrics,
                     model_selection, naive_bayes, neighbors,
                     pipeline, preprocessing, svm, linear_model,
                     tree, discriminant_analysis)
from sklearn.exceptions import FitFailedWarning
from sklearn.model_selection._validation import _score, cross_validate
from sklearn.model_selection import _search, _validation

from galaxy_ml.utils import (SafeEval, get_cv, get_scoring, load_model,
                             read_columns, try_get_attr, get_module)


_fit_and_score = try_get_attr('galaxy_ml.model_validations', '_fit_and_score')
setattr(_search, '_fit_and_score', _fit_and_score)
setattr(_validation, '_fit_and_score', _fit_and_score)

N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))
CACHE_DIR = './cached'
NON_SEARCHABLE = ('n_jobs', 'pre_dispatch', 'memory', '_path',
                  'nthread')


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


def main(inputs, infile_estimator, infile1, infile2,
         outfile_result, outfile_object=None, groups=None,
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

    params_builder = params['search_schemes']['search_params_builder']

    with open(infile_estimator, 'rb') as estimator_handler:
        estimator = load_model(estimator_handler)
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
        options['cv_selector']['groups_selector'] = groups

    splitter, groups = get_cv(options.pop('cv_selector'))
    options['cv'] = splitter
    options['n_jobs'] = N_JOBS
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

    # del loaded_df
    del loaded_df

    # handle memory
    memory = joblib.Memory(location=CACHE_DIR, verbose=0)
    # cache iraps_core fits could increase search speed significantly
    if estimator.__class__.__name__ == 'IRAPSClassifier':
        estimator.set_params(memory=memory)
    else:
        # For iraps buried in pipeline
        for p, v in estimator_params.items():
            if p.endswith('memory'):
                # for case of `__irapsclassifier__memory`
                if len(p) > 8 and p[:-8].endswith('irapsclassifier'):
                    # cache iraps_core fits could increase search
                    # speed significantly
                    new_params = {p: memory}
                    estimator.set_params(**new_params)
                # security reason, we don't want memory being
                # modified unexpectedly
                elif v:
                    new_params = {p, None}
                    estimator.set_params(**new_params)
            # For now, 1 CPU is suggested for iprasclassifier
            elif p.endswith('n_jobs'):
                new_params = {p: 1}
                estimator.set_params(**new_params)

    param_grid = _eval_search_params(params_builder)
    searcher = optimizer(estimator, param_grid, **options)

    # do nested split
    split_mode = params['outer_split'].pop('split_mode')
    # nested CV, outer cv using cross_validate
    if split_mode == 'nested_cv':
        outer_cv, _ = get_cv(params['outer_split']['cv_selector'])

        if options['error_score'] == 'raise':
            rval = cross_validate(
                searcher, X, y, scoring=options['scoring'],
                cv=outer_cv, n_jobs=N_JOBS, verbose=0,
                error_score=options['error_score'])
        else:
            warnings.simplefilter('always', FitFailedWarning)
            with warnings.catch_warnings(record=True) as w:
                try:
                    rval = cross_validate(
                        searcher, X, y,
                        scoring=options['scoring'],
                        cv=outer_cv, n_jobs=N_JOBS,
                        verbose=0,
                        error_score=options['error_score'])
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
        rval.to_csv(path_or_buf=outfile_result, sep='\t',
                    header=True, index=False)
    else:
        if split_mode == 'train_test_split':
            train_test_split = try_get_attr(
                'galaxy_ml.model_validations', 'train_test_split')
            # make sure refit is choosen
            # this could be True for sklearn models, but not the case for
            # deep learning models
            if not options['refit'] and \
                    not all(hasattr(estimator, attr)
                            for attr in ('config', 'model_type')):
                warnings.warn("Refit is change to `True` for nested "
                              "validation!")
                setattr(searcher, 'refit', True)
            split_options = params['outer_split']

            # splits
            if split_options['shuffle'] == 'stratified':
                split_options['labels'] = y
                X, X_test, y, y_test = train_test_split(X, y, **split_options)
            elif split_options['shuffle'] == 'group':
                if groups is None:
                    raise ValueError("No group based CV option was "
                                     "choosen for group shuffle!")
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
        # end train_test_split

        # shared by both train_test_split and non-split
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

        # no outer split
        if split_mode == 'no':
            # save results
            cv_results = pd.DataFrame(searcher.cv_results_)
            cv_results = cv_results[sorted(cv_results.columns)]
            cv_results.to_csv(path_or_buf=outfile_result, sep='\t',
                              header=True, index=False)

        # train_test_split, output test result using best_estimator_
        # or rebuild the trained estimator using weights if applicable.
        else:
            scorer_ = searcher.scorer_
            if isinstance(scorer_, collections.Mapping):
                is_multimetric = True
            else:
                is_multimetric = False

            best_estimator_ = getattr(searcher, 'best_estimator_', None)
            if not best_estimator_:
                weights_path = __import__('os').path.join(
                    __import__('os').getcwd(), 'weights.hdf5')
                # for keras_g deep learning models
                if __import__('os').path.exists(weights_path):
                    best_estimator_ = clone(estimator)
                    # for pipeline object
                    if isinstance(estimator, pipeline.Pipeline):
                        best_estimator_[-1][-1].load_weights(weights_path)
                    else:
                        best_estimator_.load_weights(weights_path)
                else:
                    raise ValueError("Estimator neither has `best_estimator_` "
                                     "attributes, nor outputs `weights.hdf5`!")

            if best_estimator_.__class__.__name__ == 'KerasGBatchClassifier' \
                    and hasattr(estimator.data_batch_generator, 'target_path'):
                best_estimator_.data_generator_ = \
                    best_estimator_.data_batch_generator
                best_estimator_.data_generator_.fit()
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
            result_df.to_csv(path_or_buf=outfile_result, sep='\t',
                             header=True, index=False)

    memory.clear(warn=False)

    if outfile_object:
        with open(outfile_object, 'wb') as output_handler:
            pickle.dump(searcher, output_handler, pickle.HIGHEST_PROTOCOL)


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
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile1, args.infile2,
         args.outfile_result, outfile_object=args.outfile_object,
         groups=args.groups, ref_seq=args.ref_seq, intervals=args.intervals,
         targets=args.targets, fasta_path=args.fasta_path)
