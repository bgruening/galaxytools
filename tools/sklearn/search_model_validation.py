import imblearn
import json
import numpy as np
import os
import pandas
import pickle
import skrebate
import sklearn
import sys
import xgboost
import warnings
from imblearn import under_sampling, over_sampling, combine
from imblearn.pipeline import Pipeline as imbPipeline
from sklearn import (cluster, compose, decomposition, ensemble, feature_extraction,
                    feature_selection, gaussian_process, kernel_approximation, metrics,
                    model_selection, naive_bayes, neighbors, pipeline, preprocessing,
                    svm, linear_model, tree, discriminant_analysis)
from sklearn.exceptions import FitFailedWarning
from sklearn.externals import joblib
from utils import get_cv, get_scoring, get_X_y, load_model, read_columns, SafeEval


N_JOBS = int(os.environ.get('GALAXY_SLOTS', 1))
NON_SEARCHABLE = ('n_jobs', 'pre_dispatch', 'memory', 'steps', 'nthread', 'verbose')


def get_search_params(params_builder):
    search_params = {}

    for p in params_builder['param_set']:
        search_list = p['sp_list'].strip()
        if search_list == '':
            continue

        param_name = p['sp_name']
        if param_name.lower().endswith( NON_SEARCHABLE ):
            print("Warning: `%s` is not eligible for search and was omitted!" % param_name)
            continue

        if not search_list.startswith(':'):
            safe_eval = SafeEval(load_scipy=True, load_numpy=True)
            ev = safe_eval(search_list)
            search_params[param_name] = ev
        else:
            # Have `:` before search list, asks for estimator evaluatio
            safe_eval_es = SafeEval(load_estimators=True)
            search_list = search_list[1:].strip()
            #TODO maybe add regular express check
            ev = safe_eval_es(search_list)
            preprocessors = (preprocessing.StandardScaler(), preprocessing.Binarizer(),
                            preprocessing.Imputer(), preprocessing.MaxAbsScaler(),
                            preprocessing.Normalizer(), preprocessing.MinMaxScaler(),
                            preprocessing.PolynomialFeatures(),preprocessing.RobustScaler(),
                            feature_selection.SelectKBest(), feature_selection.GenericUnivariateSelect(),
                            feature_selection.SelectPercentile(), feature_selection.SelectFpr(),
                            feature_selection.SelectFdr(), feature_selection.SelectFwe(),
                            feature_selection.VarianceThreshold(),
                            decomposition.FactorAnalysis(random_state=0),
                            decomposition.FastICA(random_state=0), decomposition.IncrementalPCA(),
                            decomposition.KernelPCA(random_state=0, n_jobs=N_JOBS),
                            decomposition.LatentDirichletAllocation(random_state=0, n_jobs=N_JOBS),
                            decomposition.MiniBatchDictionaryLearning(random_state=0, n_jobs=N_JOBS),
                            decomposition.MiniBatchSparsePCA(random_state=0, n_jobs=N_JOBS),
                            decomposition.NMF(random_state=0), decomposition.PCA(random_state=0),
                            decomposition.SparsePCA(random_state=0, n_jobs=N_JOBS),
                            decomposition.TruncatedSVD(random_state=0),
                            kernel_approximation.Nystroem(random_state=0),
                            kernel_approximation.RBFSampler(random_state=0),
                            kernel_approximation.AdditiveChi2Sampler(),
                            kernel_approximation.SkewedChi2Sampler(random_state=0),
                            cluster.FeatureAgglomeration(),
                            skrebate.ReliefF(n_jobs=N_JOBS), skrebate.SURF(n_jobs=N_JOBS),
                            skrebate.SURFstar(n_jobs=N_JOBS), skrebate.MultiSURF(n_jobs=N_JOBS),
                            skrebate.MultiSURFstar(n_jobs=N_JOBS),
                            imblearn.under_sampling.ClusterCentroids(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.CondensedNearestNeighbour(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.EditedNearestNeighbours(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.RepeatedEditedNearestNeighbours(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.AllKNN(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.InstanceHardnessThreshold(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.NearMiss(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.NeighbourhoodCleaningRule(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.OneSidedSelection(random_state=0, n_jobs=N_JOBS),
                            imblearn.under_sampling.RandomUnderSampler(random_state=0),
                            imblearn.under_sampling.TomekLinks(random_state=0, n_jobs=N_JOBS),
                            imblearn.over_sampling.ADASYN(random_state=0, n_jobs=N_JOBS),
                            imblearn.over_sampling.RandomOverSampler(random_state=0),
                            imblearn.over_sampling.SMOTE(random_state=0, n_jobs=N_JOBS),
                            imblearn.over_sampling.SVMSMOTE(random_state=0, n_jobs=N_JOBS),
                            imblearn.over_sampling.BorderlineSMOTE(random_state=0, n_jobs=N_JOBS),
                            imblearn.over_sampling.SMOTENC(categorical_features=[], random_state=0, n_jobs=N_JOBS),
                            imblearn.combine.SMOTEENN(random_state=0), imblearn.combine.SMOTETomek(random_state=0))
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
                elif  type(obj) is int and -1 < obj < len(preprocessors):
                    newlist.append(preprocessors[obj])
                elif hasattr(obj, 'get_params'):       # user uploaded object
                    if 'n_jobs' in obj.get_params():
                        newlist.append( obj.set_params(n_jobs=N_JOBS) )
                    else:
                        newlist.append(obj)
                else:
                    sys.exit("Unsupported estimator type: %r" %(obj))

            search_params[param_name] = newlist

    return search_params


if __name__ == '__main__':

    warnings.simplefilter('ignore')

    inputs = sys.argv[1].split(',')
    input_json_path = inputs[0]
    with open(input_json_path, 'r') as param_handler:
        params = json.load(param_handler)
    if len(inputs) > 1:
        params['search_schemes']['options']['cv_selector']['groups_selector']['infile_g'] = inputs[1]

    infile_estimator = sys.argv[2]
    infile1 = sys.argv[3]
    infile2 = sys.argv[4]
    outfile_result = sys.argv[5]
    if len(sys.argv) > 6:
        outfile_object = sys.argv[6]
    else:
        outfile_object = None

    params_builder = params['search_schemes']['search_params_builder']

    input_type = params['input_options']['selected_input']
    if input_type == 'tabular':
        header = 'infer' if params['input_options']['header1'] else None
        column_option = params['input_options']['column_selector_options_1']['selected_column_selector_option']
        if column_option in ['by_index_number', 'all_but_by_index_number', 'by_header_name', 'all_but_by_header_name']:
            c = params['input_options']['column_selector_options_1']['col1']
        else:
            c = None
        X = read_columns(
                infile1,
                c = c,
                c_option = column_option,
                sep='\t',
                header=header,
                parse_dates=True).astype(float)
    else:
        X = mmread(open(infile1, 'r'))

    header = 'infer' if params['input_options']['header2'] else None
    column_option = params['input_options']['column_selector_options_2']['selected_column_selector_option2']
    if column_option in ['by_index_number', 'all_but_by_index_number', 'by_header_name', 'all_but_by_header_name']:
        c = params['input_options']['column_selector_options_2']['col2']
    else:
        c = None
    y = read_columns(
            infile2,
            c = c,
            c_option = column_option,
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
        options['refit'] = 'primary'
    if 'pre_dispatch' in options and options['pre_dispatch'] == '':
        options['pre_dispatch'] = None

    with open(infile_estimator, 'rb') as estimator_handler:
        estimator = load_model(estimator_handler)

    param_grid = get_search_params(params_builder)
    searcher = optimizer(estimator, param_grid, **options)

    for p, v in estimator.get_params().items():
        if p.endswith('memory') and v:
            estimator.set_params( p=None )
        elif p.endswith('n_jobs') and v != N_JOBS:
            estimator.set_params( p=N_JOBS )

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

    cv_result = pandas.DataFrame(searcher.cv_results_)
    col_rename = {}
    for col in cv_result.columns:
        if col.endswith('_primary'):
            col_rename[col] = col[:-7] + primary_scoring
    cv_result.rename(inplace=True, columns=col_rename)
    cv_result = cv_result[sorted(cv_result.columns)]
    cv_result.to_csv(path_or_buf=outfile_result, sep='\t', header=True, index=False)

    if outfile_object:
        with open(outfile_object, 'wb') as output_handler:
            pickle.dump(searcher, output_handler, pickle.HIGHEST_PROTOCOL)
