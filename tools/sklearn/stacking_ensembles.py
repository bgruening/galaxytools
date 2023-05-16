import argparse
import ast
import json
import sys
import warnings
from distutils.version import LooseVersion as Version

from galaxy_ml import __version__ as galaxy_ml_version
from galaxy_ml.model_persist import dump_model_to_h5, load_model_from_h5
from galaxy_ml.utils import get_cv, get_estimator

import mlxtend.classifier
import mlxtend.regressor


warnings.filterwarnings('ignore')

N_JOBS = int(__import__('os').environ.get('GALAXY_SLOTS', 1))


def main(inputs_path, output_obj, base_paths=None, meta_path=None):
    """
    Parameter
    ---------
    inputs_path : str
        File path for Galaxy parameters

    output_obj : str
        File path for ensemble estimator ouput

    base_paths : str
        File path or paths concatenated by comma.

    meta_path : str
        File path
    """
    with open(inputs_path, 'r') as param_handler:
        params = json.load(param_handler)

    estimator_type = params['algo_selection']['estimator_type']
    # get base estimators
    base_estimators = []
    for idx, base_file in enumerate(base_paths.split(',')):
        if base_file and base_file != 'None':
            model = load_model_from_h5(base_file)
        else:
            estimator_json = (params['base_est_builder'][idx]
                              ['estimator_selector'])
            model = get_estimator(estimator_json)

        if estimator_type.startswith('sklearn'):
            named = model.__class__.__name__.lower()
            named = 'base_%d_%s' % (idx, named)
            base_estimators.append((named, model))
        else:
            base_estimators.append(model)

    # get meta estimator, if applicable
    if estimator_type.startswith('mlxtend'):
        if meta_path:
            meta_estimator = load_model_from_h5(meta_path)
        else:
            estimator_json = (params['algo_selection']
                              ['meta_estimator']['estimator_selector'])
            meta_estimator = get_estimator(estimator_json)

    options = params['algo_selection']['options']

    cv_selector = options.pop('cv_selector', None)
    if cv_selector:
        if Version(galaxy_ml_version) < Version('0.8.3'):
            cv_selector.pop('n_stratification_bins', None)
        splitter, groups = get_cv(cv_selector)
        options['cv'] = splitter
        # set n_jobs
        options['n_jobs'] = N_JOBS

    weights = options.pop('weights', None)
    if weights:
        weights = ast.literal_eval(weights)
        if weights:
            options['weights'] = weights

    mod_and_name = estimator_type.split('_')
    mod = sys.modules[mod_and_name[0]]
    klass = getattr(mod, mod_and_name[1])

    if estimator_type.startswith('sklearn'):
        options['n_jobs'] = N_JOBS
        ensemble_estimator = klass(base_estimators, **options)

    elif mod == mlxtend.classifier:
        ensemble_estimator = klass(
            classifiers=base_estimators,
            meta_classifier=meta_estimator,
            **options)

    else:
        ensemble_estimator = klass(
            regressors=base_estimators,
            meta_regressor=meta_estimator,
            **options)

    print(ensemble_estimator)
    for base_est in base_estimators:
        print(base_est)

    dump_model_to_h5(ensemble_estimator, output_obj)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-b", "--bases", dest="bases")
    aparser.add_argument("-m", "--meta", dest="meta")
    aparser.add_argument("-i", "--inputs", dest="inputs")
    aparser.add_argument("-o", "--outfile", dest="outfile")
    args = aparser.parse_args()

    main(args.inputs, args.outfile, base_paths=args.bases,
         meta_path=args.meta)
