import argparse
import json
import warnings
from distutils.version import LooseVersion as Version

from galaxy_ml import __version__ as galaxy_ml_version
from galaxy_ml.model_validations import train_test_split
from galaxy_ml.utils import get_cv, read_columns

import pandas as pd


def _get_single_cv_split(params, array, infile_labels=None,
                         infile_groups=None):
    """ output (train, test) subset from a cv splitter

    Parameters
    ----------
    params : dict
        Galaxy tool inputs
    array : pandas DataFrame object
        The target dataset to split
    infile_labels : str
        File path to dataset containing target values
    infile_groups : str
        File path to dataset containing group values
    """
    y = None
    groups = None

    nth_split = params['mode_selection']['nth_split']

    # read groups
    if infile_groups:
        header = 'infer' if (params['mode_selection']['cv_selector']
                             ['groups_selector']['header_g']) else None
        column_option = (params['mode_selection']['cv_selector']
                         ['groups_selector']['column_selector_options_g']
                         ['selected_column_selector_option_g'])
        if column_option in ['by_index_number', 'all_but_by_index_number',
                             'by_header_name', 'all_but_by_header_name']:
            c = (params['mode_selection']['cv_selector']['groups_selector']
                 ['column_selector_options_g']['col_g'])
        else:
            c = None

        groups = read_columns(
            infile_groups, c=c, c_option=column_option, sep='\t',
            header=header, parse_dates=True,
        )
        groups = groups.ravel()

        params['mode_selection']['cv_selector']['groups_selector'] = groups

    # read labels
    if infile_labels:
        target_input = (params['mode_selection']
                        ['cv_selector'].pop('target_input'))
        header = 'infer' if target_input['header1'] else None
        col_index = target_input['col'][0] - 1
        df = pd.read_csv(infile_labels, sep='\t', header=header,
                         parse_dates=True)
        y = df.iloc[:, col_index].values

    # construct the cv splitter object
    cv_selector = params['mode_selection']['cv_selector']
    if Version(galaxy_ml_version) < Version('0.8.3'):
        cv_selector.pop('n_stratification_bins', None)
    splitter, groups = get_cv(cv_selector)

    total_n_splits = splitter.get_n_splits(array.values, y=y, groups=groups)
    if nth_split > total_n_splits:
        raise ValueError("Total number of splits is {}, but got `nth_split` "
                         "= {}".format(total_n_splits, nth_split))

    i = 1
    for train_index, test_index in splitter.split(array.values, y=y,
                                                  groups=groups):
        # suppose nth_split >= 1
        if i == nth_split:
            break
        else:
            i += 1

    train = array.iloc[train_index, :]
    test = array.iloc[test_index, :]

    return train, test


def main(inputs, infile_array, outfile_train, outfile_test,
         infile_labels=None, infile_groups=None):
    """
    Parameter
    ---------
    inputs : str
        File path to galaxy tool parameter

    infile_array : str
        File paths of input arrays separated by comma

    infile_labels : str
        File path to dataset containing labels

    infile_groups : str
        File path to dataset containing groups

    outfile_train : str
        File path to dataset containing train split

    outfile_test : str
        File path to dataset containing test split
    """
    warnings.simplefilter('ignore')

    with open(inputs, 'r') as param_handler:
        params = json.load(param_handler)

    input_header = params['header0']
    header = 'infer' if input_header else None
    array = pd.read_csv(infile_array, sep='\t', header=header,
                        parse_dates=True)

    # train test split
    if params['mode_selection']['selected_mode'] == 'train_test_split':
        options = params['mode_selection']['options']
        shuffle_selection = options.pop('shuffle_selection')
        options['shuffle'] = shuffle_selection['shuffle']
        if infile_labels:
            header = 'infer' if shuffle_selection['header1'] else None
            col_index = shuffle_selection['col'][0] - 1
            df = pd.read_csv(infile_labels, sep='\t', header=header,
                             parse_dates=True)
            labels = df.iloc[:, col_index].values
            options['labels'] = labels

        train, test = train_test_split(array, **options)

    # cv splitter
    else:
        train, test = _get_single_cv_split(params, array,
                                           infile_labels=infile_labels,
                                           infile_groups=infile_groups)

    print("Input shape: %s" % repr(array.shape))
    print("Train shape: %s" % repr(train.shape))
    print("Test shape: %s" % repr(test.shape))
    train.to_csv(outfile_train, sep='\t', header=input_header, index=False)
    test.to_csv(outfile_test, sep='\t', header=input_header, index=False)


if __name__ == '__main__':
    aparser = argparse.ArgumentParser()
    aparser.add_argument("-i", "--inputs", dest="inputs", required=True)
    aparser.add_argument("-X", "--infile_array", dest="infile_array")
    aparser.add_argument("-y", "--infile_labels", dest="infile_labels")
    aparser.add_argument("-g", "--infile_groups", dest="infile_groups")
    aparser.add_argument("-o", "--outfile_train", dest="outfile_train")
    aparser.add_argument("-t", "--outfile_test", dest="outfile_test")
    args = aparser.parse_args()

    main(args.inputs, args.infile_array, args.outfile_train,
         args.outfile_test, args.infile_labels, args.infile_groups)
