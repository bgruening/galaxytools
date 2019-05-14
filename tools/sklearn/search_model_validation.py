import argparse
import collections
import imblearn
import json
import numpy as np
import pandas
import pickle
import skrebate
import sklearn
import sys
import xgboost
import warnings
import iraps_classifier
import model_validations
import preprocessors
import feature_selectors
from imblearn import under_sampling, over_sampling, combine
from scipy.io import mmread
from mlxtend import classifier, regressor
from sklearn import (cluster, compose, decomposition, ensemble,
                     feature_extraction, feature_selection,
                     gaussian_process, kernel_approximation, metrics,
                     model_selection, naive_bayes, neighbors,
                     pipeline, preprocessing, svm, linear_model,
                     tree, discriminant_analysis)
from sklearn.exceptions import FitFailedWarning
from sklearn.externals import joblib
from sklearn.model_selection._validation import _score

from utils import (SafeEval, get_cv, get_scoring, get_X_y,
                   load_model, read_columns)
from model_validations import train_test_split


N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))
CACHE_DIR = './cached'
NON_SEARCHABLE = ('n_jobs', 'pre_dispatch', 'memory', 'steps',
                  'nthread', 'verbose')


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
            preprocessors = (
                preprocessing.StandardScaler(), preprocessing.Binarizer(),
                preprocessing.Imputer(), preprocessing.MaxAbsScaler(),
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
                    newlist.extend(preprocessors[0:36])
                elif obj == 'sk_prep_all':      # no KernalCenter()
                    newlist.extend(preprocessors[0:8])
                elif obj == 'fs_all':
                    newlist.extend(preprocessors[8:15])
                elif obj == 'decomp_all':
                    newlist.extend(preprocessors[15:26])
                elif obj == 'k_appr_all':
                    newlist.extend(preprocessors[26:30])
                elif obj == 'reb_all':
                    newlist.extend(preprocessors[31:36])
                elif obj == 'imb_all':
                    newlist.extend(preprocessors[36:55])
                elif type(obj) is int and -1 < obj < len(preprocessors):
                    newlist.append(preprocessors[obj])
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
         outfile_result, outfile_object=None, groups=None):
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
    """

    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)
    if groups:
        (params['search_schemes']['options']['cv_selector']
         ['groups_selector']['infile_g']) = groups

    params_builder = params['search_schemes']['search_params_builder']

    input_type = params['input_options']['selected_input']
    if input_type == 'tabular':
        header = 'infer' if params['input_options']['header1'] else None
        column_option = (params['input_options']['column_selector_options_1']
                         ['selected_column_selector_option'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = params['input_options']['column_selector_options_1']['col1']
        else:
            c = None
        X = read_columns(
                infile1,
                c=c,
                c_option=column_option,
                sep='\t',
                header=header,
                parse_dates=True).astype(float)
    else:
        X = mmread(open(infile1, 'r'))

    header = 'infer' if params['input_options']['header2'] else None
    column_option = (params['input_options']['column_selector_options_2']
                     ['selected_column_selector_option2'])
    if column_option in ['by_index_number', 'all_but_by_index_number',
                         'by_header_name', 'all_but_by_header_name']:
        c = params['input_options']['column_selector_options_2']['col2']
    else:
        c = None
    y = read_columns(
            infile2,
            c=c,
            c_option=column_option,
            sep='\t',
            header=header,
            parse_dates=True)
    y = y.ravel()

    optimizer = params['search_schemes']['selected_search_scheme']
    optimizer = getattr(model_selection, optimizer)

    options = params['search_schemes']['options']

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

    with open(infile_estimator, 'rb') as estimator_handler:
        estimator = load_model(estimator_handler)

    memory = joblib.Memory(location=CACHE_DIR, verbose=0)
    # cache iraps_core fits could increase search speed significantly
    if estimator.__class__.__name__ == 'IRAPSClassifier':
        estimator.set_params(memory=memory)
    else:
        for p, v in estimator.get_params().items():
            if p.endswith('memory'):
                if len(p) > 8 and p[:-8].endswith('irapsclassifier'):
                    # cache iraps_core fits could increase search
                    # speed significantly
                    new_params = {p: memory}
                    estimator.set_params(**new_params)
                elif v:
                    new_params = {p, None}
                    estimator.set_params(**new_params)
            elif p.endswith('n_jobs'):
                new_params = {p: 1}
                estimator.set_params(**new_params)

    param_grid = _eval_search_params(params_builder)
    searcher = optimizer(estimator, param_grid, **options)

    # do train_test_split
    do_train_test_split = params['train_test_split'].pop('do_split')
    if do_train_test_split == 'yes':
        # make sure refit is choosen
        if not options['refit']:
            raise ValueError("Refit must be `True` for shuffle splitting!")
        split_options = params['train_test_split']

        # splits
        if split_options['shuffle'] == 'stratified':
            split_options['labels'] = y
            X, X_test, y, y_test = train_test_split(X, y, **split_options)
        elif split_options['shuffle'] == 'group':
            if not groups:
                raise ValueError("No group based CV option was "
                                 "choosen for group shuffle!")
            split_options['labels'] = groups
            X, X_test, y, y_test, groups, _ =\
                train_test_split(X, y, **split_options)
        else:
            if split_options['shuffle'] == 'None':
                split_options['shuffle'] = None
            X, X_test, y, y_test =\
                train_test_split(X, y, **split_options)
    # end train_test_split

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

    if do_train_test_split == 'no':
        # save results
        cv_results = pandas.DataFrame(searcher.cv_results_)
        cv_results = cv_results[sorted(cv_results.columns)]
        cv_results.to_csv(path_or_buf=outfile_result, sep='\t',
                          header=True, index=False)

    # output test result using best_estimator_
    else:
        best_estimator_ = searcher.best_estimator_
        if isinstance(options['scoring'], collections.Mapping):
            is_multimetric = True
        else:
            is_multimetric = False

        test_score = _score(best_estimator_, X_test,
                            y_test, options['scoring'],
                            is_multimetric=is_multimetric)
        if not is_multimetric:
            test_score = {primary_scoring: test_score}
        for key, value in test_score.items():
            test_score[key] = [value]
        result_df = pandas.DataFrame(test_score)
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
    aparser.add_argument("-r", "--outfile_result", dest="outfile_result")
    aparser.add_argument("-o", "--outfile_object", dest="outfile_object")
    aparser.add_argument("-g", "--groups", dest="groups")
    args = aparser.parse_args()

    main(args.inputs, args.infile_estimator, args.infile1, args.infile2,
         args.outfile_result, outfile_object=args.outfile_object,
         groups=args.groups)
