import sys
import os
import pandas
import re
import pickle
import numpy as np
import xgboost
import scipy
import sklearn
import ast
from asteval import Interpreter, make_symbol_table
from sklearn import metrics, model_selection, ensemble, svm, linear_model, naive_bayes, tree, neighbors


def read_columns(f, c=None, c_option='by_index_number', return_df=False, **args):
  data = pandas.read_csv(f, **args)
  if c_option == 'by_index_number':
    cols = list(map(lambda x: x - 1, c))
    data = data.iloc[:,cols]
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
  return y


## generate an instance for one of sklearn.feature_selection classes
def feature_selector(inputs):
  selector = inputs["selected_algorithm"]
  selector = getattr(sklearn.feature_selection, selector)
  options = inputs["options"]

  if inputs['selected_algorithm'] == 'SelectFromModel':
    if not options['threshold'] or options['threshold'] == 'None':
      options['threshold'] = None
    if inputs['model_inputter']['input_mode'] == 'prefitted':
      model_file = inputs['model_inputter']['fitted_estimator']
      with open(model_file, 'rb') as model_handler:
        fitted_estimator = pickle.load(model_handler)
      new_selector = selector(fitted_estimator, prefit=True, **options)
    else:
      estimator_json = inputs['model_inputter']["estimator_selector"]
      estimator = get_estimator(estimator_json)
      new_selector = selector(estimator, **options)

  elif inputs['selected_algorithm'] == 'RFE':
    estimator=get_estimator(inputs["estimator_selector"])
    new_selector = selector(estimator, **options)

  elif inputs['selected_algorithm'] == 'RFECV':
    options['scoring'] = get_scoring(options['scoring'])
    options['n_jobs'] = int( os.environ.get('GALAXY_SLOTS', 1) )
    options['cv'] = get_cv( options['cv'].strip() )
    estimator=get_estimator(inputs["estimator_selector"])
    new_selector = selector(estimator, **options)

  elif inputs['selected_algorithm'] == "VarianceThreshold":
    new_selector = selector(**options)

  else:
    score_func = inputs["score_func"]
    score_func = getattr(sklearn.feature_selection, score_func)
    new_selector = selector(score_func, **options)

  return new_selector
 

def get_X_y(params, file1, file2):
  input_type = params["selected_tasks"]["selected_algorithms"]["input_options"]["selected_input"]
  if input_type=="tabular":
    header = 'infer' if params["selected_tasks"]["selected_algorithms"]["input_options"]["header1"] else None
    column_option = params["selected_tasks"]["selected_algorithms"]["input_options"]["column_selector_options_1"]["selected_column_selector_option"]
    if column_option in ["by_index_number", "all_but_by_index_number", "by_header_name", "all_but_by_header_name"]:
      c = params["selected_tasks"]["selected_algorithms"]["input_options"]["column_selector_options_1"]["col1"]
    else:
      c = None
    X = read_columns(
      file1,
      c = c,
      c_option = column_option,
      sep='\t',
      header=header,
      parse_dates=True
    )
  else:
    X = mmread(file1)

  header = 'infer' if params["selected_tasks"]["selected_algorithms"]["input_options"]["header2"] else None
  column_option = params["selected_tasks"]["selected_algorithms"]["input_options"]["column_selector_options_2"]["selected_column_selector_option2"]
  if column_option in ["by_index_number", "all_but_by_index_number", "by_header_name", "all_but_by_header_name"]:
    c = params["selected_tasks"]["selected_algorithms"]["input_options"]["column_selector_options_2"]["col2"]
  else:
    c = None
  y = read_columns(
    file2,
    c = c,
    c_option = column_option,
    sep='\t',
    header=header,
    parse_dates=True
  )
  y=y.ravel()
  return X, y


class SafeEval(Interpreter):

  def __init__(self, load_scipy=False, load_numpy=False):

    # File opening and other unneeded functions could be dropped
    unwanted = ['open', 'type', 'dir', 'id', 'str', 'repr']

    # Allowed symbol table. Add more if needed.
    new_syms = {
      'np_arange': getattr(np, 'arange'),
      'ensemble_ExtraTreesClassifier': getattr(ensemble, 'ExtraTreesClassifier')
    }

    syms = make_symbol_table(use_numpy=False, **new_syms)

    if load_scipy:
      scipy_distributions = scipy.stats.distributions.__dict__
      for key in scipy_distributions.keys():
        if isinstance(scipy_distributions[key], (scipy.stats.rv_continuous, scipy.stats.rv_discrete)):
          syms['scipy_stats_' + key] = scipy_distributions[key]

    if load_numpy:
      from_numpy_random = [ 'beta', 'binomial', 'bytes', 'chisquare', 'choice', 'dirichlet', 'division',
                            'exponential', 'f', 'gamma', 'geometric', 'gumbel', 'hypergeometric',
                            'laplace', 'logistic', 'lognormal', 'logseries', 'mtrand', 'multinomial',
                            'multivariate_normal', 'negative_binomial', 'noncentral_chisquare', 'noncentral_f',
                            'normal', 'pareto', 'permutation', 'poisson', 'power', 'rand', 'randint',
                            'randn', 'random', 'random_integers', 'random_sample', 'ranf', 'rayleigh',
                            'sample', 'seed', 'set_state', 'shuffle', 'standard_cauchy', 'standard_exponential',
                            'standard_gamma', 'standard_normal', 'standard_t', 'triangular', 'uniform',
                            'vonmises', 'wald', 'weibull', 'zipf' ]
      for f in from_numpy_random:
        syms['np_random_' + f] = getattr(np.random, f)

    for key in unwanted:
      syms.pop(key, None)

    super(SafeEval, self).__init__(symtable=syms, use_numpy=False, minimal=False,
                                  no_if=True, no_for=True, no_while=True, no_try=True,
                                  no_functiondef=True, no_ifexp=True, no_listcomp=False,
                                  no_augassign=False, no_assert=True, no_delete=True,
                                  no_raise=True, no_print=True)


def get_search_params(params_builder):
  search_params = {}
  safe_eval = SafeEval(load_scipy=True, load_numpy=True)

  for p in params_builder['param_set']:
    search_p = p['search_param_selector']['search_p']
    if search_p.strip() == '':
      continue
    param_type = p['search_param_selector']['selected_param_type']

    lst = search_p.split(":")
    assert (len(lst) == 2), "Error, make sure there is one and only one colon in search parameter input."
    literal = lst[1].strip()
    ev = safe_eval(literal)
    if param_type == "final_estimator_p":
      search_params["estimator__" + lst[0].strip()] = ev
    else:
      search_params["preprocessing_" + param_type[5:6] + "__" + lst[0].strip()] = ev

  return search_params


def get_estimator(estimator_json):
  estimator_module = estimator_json['selected_module']
  estimator_cls = estimator_json['selected_estimator']

  if estimator_module == "xgboost":
    cls = getattr(xgboost, estimator_cls)
  else:
    module = getattr(sklearn, estimator_module)
    cls = getattr(module, estimator_cls)

  estimator = cls()

  estimator_params = estimator_json['text_params'].strip()
  if estimator_params != "":
    try:
      params = ast.literal_eval('{' + estimator_params + '}')
    except ValueError:
      sys.exit("Unsupported parameter input: `%s`" %estimator_params)
    estimator.set_params(**params)

  return estimator


def get_cv(literal):
  safe_eval = SafeEval()
  if literal == "":
    return None
  if literal.isdigit():
    return int(literal)
  m = re.match(r'^(?P<method>\w+)\((?P<args>.*)\)$', literal)
  if m:
    my_class = getattr( model_selection, m.group('method') )
    args = safe_eval( 'dict('+ m.group('args') + ')' )
    return my_class( **args )
  sys.exit("Unsupported CV input: %s" %literal)


def get_scoring(scorer):
  import warnings
  def balanced_accuracy_score(y_true, y_pred):
      C = metrics.confusion_matrix(y_true, y_pred)
      with np.errstate(divide='ignore', invalid='ignore'):
          per_class = np.diag(C) / C.sum(axis=1)
      if np.any(np.isnan(per_class)):
          warnings.warn('y_pred contains classes not in y_true')
          per_class = per_class[~np.isnan(per_class)]
      score = np.mean(per_class)
      return score

  if scorer == "default":
    return None
  if scorer == "balanced_accuracy":
    sk_version = tuple([int(x) for x in sklearn.__version__.split('.')])
    if sk_version < (0, 20):
      return metrics.make_scorer(balanced_accuracy_score)
    else:
      return scorer
  return scorer

