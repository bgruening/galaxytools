import ast
import json
import imblearn
import numpy as np
import pandas
import pickle
import re
import scipy
import sklearn
import skrebate
import sys
import warnings
import xgboost

from collections import Counter
from asteval import Interpreter, make_symbol_table
from imblearn import under_sampling, over_sampling, combine
from imblearn.pipeline import Pipeline as imbPipeline
from mlxtend import regressor, classifier
from scipy.io import mmread
from sklearn import (
    cluster, compose, decomposition, ensemble, feature_extraction,
    feature_selection, gaussian_process, kernel_approximation, metrics,
    model_selection, naive_bayes, neighbors, pipeline, preprocessing,
    svm, linear_model, tree, discriminant_analysis)

try:
    import iraps_classifier
except ImportError:
    pass

try:
    import model_validations
except ImportError:
    pass

try:
    import feature_selectors
except ImportError:
    pass

try:
    import preprocessors
except ImportError:
    pass

# handle pickle white list file
WL_FILE = __import__('os').path.join(
    __import__('os').path.dirname(__file__), 'pk_whitelist.json')

N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))


class _SafePickler(pickle.Unpickler, object):
    """
    Used to safely deserialize scikit-learn model objects
    Usage:
        eg.: _SafePickler.load(pickled_file_object)
    """
    def __init__(self, file):
        super(_SafePickler, self).__init__(file)
        # load global white list
        with open(WL_FILE, 'r') as f:
            self.pk_whitelist = json.load(f)

        self.bad_names = (
            'and', 'as', 'assert', 'break', 'class', 'continue',
            'def', 'del', 'elif', 'else', 'except', 'exec',
            'finally', 'for', 'from', 'global', 'if', 'import',
            'in', 'is', 'lambda', 'not', 'or', 'pass', 'print',
            'raise', 'return', 'try', 'system', 'while', 'with',
            'True', 'False', 'None', 'eval', 'execfile', '__import__',
            '__package__', '__subclasses__', '__bases__', '__globals__',
            '__code__', '__closure__', '__func__', '__self__', '__module__',
            '__dict__', '__class__', '__call__', '__get__',
            '__getattribute__', '__subclasshook__', '__new__',
            '__init__', 'func_globals', 'func_code', 'func_closure',
            'im_class', 'im_func', 'im_self', 'gi_code', 'gi_frame',
            '__asteval__', 'f_locals', '__mro__')

        # unclassified good globals
        self.good_names = [
            'copy_reg._reconstructor', '__builtin__.object',
            '__builtin__.bytearray', 'builtins.object',
            'builtins.bytearray', 'keras.engine.sequential.Sequential',
            'keras.engine.sequential.Model']

        # custom module in Galaxy-ML
        self.custom_modules = [
            '__main__', 'keras_galaxy_models', 'feature_selectors',
            'preprocessors', 'iraps_classifier', 'model_validations']

    # override
    def find_class(self, module, name):
        # balack list first
        if name in self.bad_names:
            raise pickle.UnpicklingError("global '%s.%s' is forbidden"
                                         % (module, name))

        # custom module in Galaxy-ML
        if module in self.custom_modules:
            cutom_module = sys.modules.get(module, None)
            if cutom_module:
                return getattr(cutom_module, name)
            else:
                raise pickle.UnpicklingError("Module %s' is not imported"
                                             % module)

        # For objects from outside libraries, it's necessary to verify
        # both module and name. Currently only a blacklist checker
        # is working.
        # TODO: replace with a whitelist checker.
        good_names = self.good_names
        pk_whitelist = self.pk_whitelist
        if re.match(r'^[a-zA-Z_][a-zA-Z0-9_]*$', name):
            fullname = module + '.' + name
            if (fullname in good_names)\
                or (module.startswith(('sklearn.', 'xgboost.', 'skrebate.',
                                       'imblearn.', 'mlxtend.', 'numpy.'))
                    or module == 'numpy'):
                if fullname not in (pk_whitelist['SK_NAMES'] +
                                    pk_whitelist['SKR_NAMES'] +
                                    pk_whitelist['XGB_NAMES'] +
                                    pk_whitelist['NUMPY_NAMES'] +
                                    pk_whitelist['IMBLEARN_NAMES'] +
                                    pk_whitelist['MLXTEND_NAMES'] +
                                    good_names):
                    # raise pickle.UnpicklingError
                    print("Warning: global %s is not in pickler whitelist "
                          "yet and will loss support soon. Contact tool "
                          "author or leave a message at github.com" % fullname)
                mod = sys.modules[module]
                return getattr(mod, name)

        raise pickle.UnpicklingError("global '%s' is forbidden" % fullname)


def load_model(file):
    """Load pickled object with `_SafePicker`
    """
    return _SafePickler(file).load()


def read_columns(f, c=None, c_option='by_index_number',
                 return_df=False, **args):
    """Return array from a tabular dataset by various columns selection
    """
    data = pandas.read_csv(f, **args)
    if c_option == 'by_index_number':
        cols = list(map(lambda x: x - 1, c))
        data = data.iloc[:, cols]
    if c_option == 'all_but_by_index_number':
        cols = list(map(lambda x: x - 1, c))
        data.drop(data.columns[cols], axis=1, inplace=True)
    if c_option == 'by_header_name':
        cols = [e.strip() for e in c.split(',')]
        data = data[cols]
    if c_option == 'all_but_by_header_name':
        cols = [e.strip() for e in c.split(',')]
        data.drop(cols, axis=1, inplace=True)
    y = data.values
    if return_df:
        return y, data
    else:
        return y


def feature_selector(inputs, X=None, y=None):
    """generate an instance of sklearn.feature_selection classes

    Parameters
    ----------
    inputs : dict
        From galaxy tool parameters.
    X : array
        Containing training features.
    y : array or list
        Target values.
    """
    selector = inputs['selected_algorithm']
    if selector != 'DyRFECV':
        selector = getattr(sklearn.feature_selection, selector)
    options = inputs['options']

    if inputs['selected_algorithm'] == 'SelectFromModel':
        if not options['threshold'] or options['threshold'] == 'None':
            options['threshold'] = None
        else:
            try:
                options['threshold'] = float(options['threshold'])
            except ValueError:
                pass
        if inputs['model_inputter']['input_mode'] == 'prefitted':
            model_file = inputs['model_inputter']['fitted_estimator']
            with open(model_file, 'rb') as model_handler:
                fitted_estimator = load_model(model_handler)
            new_selector = selector(fitted_estimator, prefit=True, **options)
        else:
            estimator_json = inputs['model_inputter']['estimator_selector']
            estimator = get_estimator(estimator_json)
            check_feature_importances = try_get_attr(
                'feature_selectors', 'check_feature_importances')
            estimator = check_feature_importances(estimator)
            new_selector = selector(estimator, **options)

    elif inputs['selected_algorithm'] == 'RFE':
        step = options.get('step', None)
        if step and step >= 1.0:
            options['step'] = int(step)
        estimator = get_estimator(inputs["estimator_selector"])
        check_feature_importances = try_get_attr(
            'feature_selectors', 'check_feature_importances')
        estimator = check_feature_importances(estimator)
        new_selector = selector(estimator, **options)

    elif inputs['selected_algorithm'] == 'RFECV':
        options['scoring'] = get_scoring(options['scoring'])
        options['n_jobs'] = N_JOBS
        splitter, groups = get_cv(options.pop('cv_selector'))
        if groups is None:
            options['cv'] = splitter
        else:
            options['cv'] = list(splitter.split(X, y, groups=groups))
        step = options.get('step', None)
        if step and step >= 1.0:
            options['step'] = int(step)
        estimator = get_estimator(inputs['estimator_selector'])
        check_feature_importances = try_get_attr(
            'feature_selectors', 'check_feature_importances')
        estimator = check_feature_importances(estimator)
        new_selector = selector(estimator, **options)

    elif inputs['selected_algorithm'] == 'DyRFECV':
        options['scoring'] = get_scoring(options['scoring'])
        options['n_jobs'] = N_JOBS
        splitter, groups = get_cv(options.pop('cv_selector'))
        if groups is None:
            options['cv'] = splitter
        else:
            options['cv'] = list(splitter.split(X, y, groups=groups))
        step = options.get('step')
        if not step or step == 'None':
            step = None
        else:
            step = ast.literal_eval(step)
        options['step'] = step
        estimator = get_estimator(inputs["estimator_selector"])
        check_feature_importances = try_get_attr(
            'feature_selectors', 'check_feature_importances')
        estimator = check_feature_importances(estimator)
        DyRFECV = try_get_attr('feature_selectors', 'DyRFECV')

        new_selector = DyRFECV(estimator, **options)

    elif inputs['selected_algorithm'] == 'VarianceThreshold':
        new_selector = selector(**options)

    else:
        score_func = inputs['score_func']
        score_func = getattr(sklearn.feature_selection, score_func)
        new_selector = selector(score_func, **options)

    return new_selector


def get_X_y(params, file1, file2):
    """Return machine learning inputs X, y from tabluar inputs
    """
    input_type = (params['selected_tasks']['selected_algorithms']
                  ['input_options']['selected_input'])
    if input_type == 'tabular':
        header = 'infer' if (params['selected_tasks']['selected_algorithms']
                             ['input_options']['header1']) else None
        column_option = (params['selected_tasks']['selected_algorithms']
                         ['input_options']['column_selector_options_1']
                         ['selected_column_selector_option'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = (params['selected_tasks']['selected_algorithms']
                 ['input_options']['column_selector_options_1']['col1'])
        else:
            c = None
        X = read_columns(
            file1,
            c=c,
            c_option=column_option,
            sep='\t',
            header=header,
            parse_dates=True).astype(float)
    else:
        X = mmread(file1)

    header = 'infer' if (params['selected_tasks']['selected_algorithms']
                         ['input_options']['header2']) else None
    column_option = (params['selected_tasks']['selected_algorithms']
                     ['input_options']['column_selector_options_2']
                     ['selected_column_selector_option2'])
    if column_option in ['by_index_number', 'all_but_by_index_number',
                         'by_header_name', 'all_but_by_header_name']:
        c = (params['selected_tasks']['selected_algorithms']
             ['input_options']['column_selector_options_2']['col2'])
    else:
        c = None
    y = read_columns(
        file2,
        c=c,
        c_option=column_option,
        sep='\t',
        header=header,
        parse_dates=True)
    y = y.ravel()

    return X, y


class SafeEval(Interpreter):
    """Customized symbol table for safely literal eval
    """
    def __init__(self, load_scipy=False, load_numpy=False,
                 load_estimators=False):

        # File opening and other unneeded functions could be dropped
        unwanted = ['open', 'type', 'dir', 'id', 'str', 'repr']

        # Allowed symbol table. Add more if needed.
        new_syms = {
            'np_arange': getattr(np, 'arange'),
            'ensemble_ExtraTreesClassifier':
                getattr(ensemble, 'ExtraTreesClassifier')
        }

        syms = make_symbol_table(use_numpy=False, **new_syms)

        if load_scipy:
            scipy_distributions = scipy.stats.distributions.__dict__
            for k, v in scipy_distributions.items():
                if isinstance(v, (scipy.stats.rv_continuous,
                                  scipy.stats.rv_discrete)):
                    syms['scipy_stats_' + k] = v

        if load_numpy:
            from_numpy_random = [
                'beta', 'binomial', 'bytes', 'chisquare', 'choice',
                'dirichlet', 'division', 'exponential', 'f', 'gamma',
                'geometric', 'gumbel', 'hypergeometric', 'laplace',
                'logistic', 'lognormal', 'logseries', 'mtrand',
                'multinomial', 'multivariate_normal', 'negative_binomial',
                'noncentral_chisquare', 'noncentral_f', 'normal', 'pareto',
                'permutation', 'poisson', 'power', 'rand', 'randint',
                'randn', 'random', 'random_integers', 'random_sample',
                'ranf', 'rayleigh', 'sample', 'seed', 'set_state',
                'shuffle', 'standard_cauchy', 'standard_exponential',
                'standard_gamma', 'standard_normal', 'standard_t',
                'triangular', 'uniform', 'vonmises', 'wald', 'weibull', 'zipf']
            for f in from_numpy_random:
                syms['np_random_' + f] = getattr(np.random, f)

        if load_estimators:
            estimator_table = {
                'sklearn_svm': getattr(sklearn, 'svm'),
                'sklearn_tree': getattr(sklearn, 'tree'),
                'sklearn_ensemble': getattr(sklearn, 'ensemble'),
                'sklearn_neighbors': getattr(sklearn, 'neighbors'),
                'sklearn_naive_bayes': getattr(sklearn, 'naive_bayes'),
                'sklearn_linear_model': getattr(sklearn, 'linear_model'),
                'sklearn_cluster': getattr(sklearn, 'cluster'),
                'sklearn_decomposition': getattr(sklearn, 'decomposition'),
                'sklearn_preprocessing': getattr(sklearn, 'preprocessing'),
                'sklearn_feature_selection':
                    getattr(sklearn, 'feature_selection'),
                'sklearn_kernel_approximation':
                    getattr(sklearn, 'kernel_approximation'),
                'skrebate_ReliefF': getattr(skrebate, 'ReliefF'),
                'skrebate_SURF': getattr(skrebate, 'SURF'),
                'skrebate_SURFstar': getattr(skrebate, 'SURFstar'),
                'skrebate_MultiSURF': getattr(skrebate, 'MultiSURF'),
                'skrebate_MultiSURFstar': getattr(skrebate, 'MultiSURFstar'),
                'skrebate_TuRF': getattr(skrebate, 'TuRF'),
                'xgboost_XGBClassifier': getattr(xgboost, 'XGBClassifier'),
                'xgboost_XGBRegressor': getattr(xgboost, 'XGBRegressor'),
                'imblearn_over_sampling': getattr(imblearn, 'over_sampling'),
                'imblearn_combine': getattr(imblearn, 'combine')
            }
            syms.update(estimator_table)

        for key in unwanted:
            syms.pop(key, None)

        super(SafeEval, self).__init__(
            symtable=syms, use_numpy=False, minimal=False,
            no_if=True, no_for=True, no_while=True, no_try=True,
            no_functiondef=True, no_ifexp=True, no_listcomp=False,
            no_augassign=False, no_assert=True, no_delete=True,
            no_raise=True, no_print=True)


def get_estimator(estimator_json):
    """Return a sklearn or compatible estimator from Galaxy tool inputs
    """
    estimator_module = estimator_json['selected_module']

    if estimator_module == 'custom_estimator':
        c_estimator = estimator_json['c_estimator']
        with open(c_estimator, 'rb') as model_handler:
            new_model = load_model(model_handler)
        return new_model

    if estimator_module == "binarize_target":
        wrapped_estimator = estimator_json['wrapped_estimator']
        with open(wrapped_estimator, 'rb') as model_handler:
            wrapped_estimator = load_model(model_handler)
        options = {}
        if estimator_json['z_score'] is not None:
            options['z_score'] = estimator_json['z_score']
        if estimator_json['value'] is not None:
            options['value'] = estimator_json['value']
        options['less_is_positive'] = estimator_json['less_is_positive']
        if estimator_json['clf_or_regr'] == 'BinarizeTargetClassifier':
            klass = try_get_attr('iraps_classifier',
                                 'BinarizeTargetClassifier')
        else:
            klass = try_get_attr('iraps_classifier',
                                 'BinarizeTargetRegressor')
        return klass(wrapped_estimator, **options)

    estimator_cls = estimator_json['selected_estimator']

    if estimator_module == 'xgboost':
        klass = getattr(xgboost, estimator_cls)
    else:
        module = getattr(sklearn, estimator_module)
        klass = getattr(module, estimator_cls)

    estimator = klass()

    estimator_params = estimator_json['text_params'].strip()
    if estimator_params != '':
        try:
            safe_eval = SafeEval()
            params = safe_eval('dict(' + estimator_params + ')')
        except ValueError:
            sys.exit("Unsupported parameter input: `%s`" % estimator_params)
        estimator.set_params(**params)
    if 'n_jobs' in estimator.get_params():
        estimator.set_params(n_jobs=N_JOBS)

    return estimator


def get_cv(cv_json):
    """ Return CV splitter from Galaxy tool inputs

    Parameters
    ----------
    cv_json : dict
        From Galaxy tool inputs.
        e.g.:
            {
                'selected_cv': 'StratifiedKFold',
                'n_splits': 3,
                'shuffle': True,
                'random_state': 0
            }
    """
    cv = cv_json.pop('selected_cv')
    if cv == 'default':
        return cv_json['n_splits'], None

    groups = cv_json.pop('groups_selector', None)
    if groups is not None:
        infile_g = groups['infile_g']
        header = 'infer' if groups['header_g'] else None
        column_option = (groups['column_selector_options_g']
                         ['selected_column_selector_option_g'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = groups['column_selector_options_g']['col_g']
        else:
            c = None
        groups = read_columns(
                infile_g,
                c=c,
                c_option=column_option,
                sep='\t',
                header=header,
                parse_dates=True)
        groups = groups.ravel()

    for k, v in cv_json.items():
        if v == '':
            cv_json[k] = None

    test_fold = cv_json.get('test_fold', None)
    if test_fold:
        if test_fold.startswith('__ob__'):
            test_fold = test_fold[6:]
        if test_fold.endswith('__cb__'):
            test_fold = test_fold[:-6]
        cv_json['test_fold'] = [int(x.strip()) for x in test_fold.split(',')]

    test_size = cv_json.get('test_size', None)
    if test_size and test_size > 1.0:
        cv_json['test_size'] = int(test_size)

    if cv == 'OrderedKFold':
        cv_class = try_get_attr('model_validations', 'OrderedKFold')
    elif cv == 'RepeatedOrderedKFold':
        cv_class = try_get_attr('model_validations', 'RepeatedOrderedKFold')
    else:
        cv_class = getattr(model_selection, cv)
    splitter = cv_class(**cv_json)

    return splitter, groups


# needed when sklearn < v0.20
def balanced_accuracy_score(y_true, y_pred):
    """Compute balanced accuracy score, which is now available in
        scikit-learn from v0.20.0.
    """
    C = metrics.confusion_matrix(y_true, y_pred)
    with np.errstate(divide='ignore', invalid='ignore'):
        per_class = np.diag(C) / C.sum(axis=1)
    if np.any(np.isnan(per_class)):
        warnings.warn('y_pred contains classes not in y_true')
        per_class = per_class[~np.isnan(per_class)]
    score = np.mean(per_class)
    return score


def get_scoring(scoring_json):
    """Return single sklearn scorer class
        or multiple scoers in dictionary
    """
    if scoring_json['primary_scoring'] == 'default':
        return None

    my_scorers = metrics.SCORERS
    my_scorers['binarize_auc_scorer'] =\
        try_get_attr('iraps_classifier', 'binarize_auc_scorer')
    my_scorers['binarize_average_precision_scorer'] =\
        try_get_attr('iraps_classifier', 'binarize_average_precision_scorer')
    if 'balanced_accuracy' not in my_scorers:
        my_scorers['balanced_accuracy'] =\
            metrics.make_scorer(balanced_accuracy_score)

    if scoring_json['secondary_scoring'] != 'None'\
            and scoring_json['secondary_scoring'] !=\
            scoring_json['primary_scoring']:
        return_scoring = {}
        primary_scoring = scoring_json['primary_scoring']
        return_scoring[primary_scoring] = my_scorers[primary_scoring]
        for scorer in scoring_json['secondary_scoring'].split(','):
            if scorer != scoring_json['primary_scoring']:
                return_scoring[scorer] = my_scorers[scorer]
        return return_scoring

    return my_scorers[scoring_json['primary_scoring']]


def get_search_params(estimator):
    """Format the output of `estimator.get_params()`
    """
    params = estimator.get_params()
    results = []
    for k, v in params.items():
        # params below won't be shown for search in the searchcv tool
        keywords = ('n_jobs', 'pre_dispatch', 'memory', 'steps',
                    'nthread', 'verbose')
        if k.endswith(keywords):
            results.append(['*', k, k+": "+repr(v)])
        else:
            results.append(['@', k, k+": "+repr(v)])
    results.append(
        ["", "Note:",
         "@, params eligible for search in searchcv tool."])

    return results


def try_get_attr(module, name):
    """try to get attribute from a custom module

    Parameters
    ----------
    module : str
        Module name
    name : str
        Attribute (class/function) name.

    Returns
    -------
    class or function
    """
    mod = sys.modules.get(module, None)
    if mod:
        return getattr(mod, name)
    else:
        raise Exception("No module named %s." % module)
