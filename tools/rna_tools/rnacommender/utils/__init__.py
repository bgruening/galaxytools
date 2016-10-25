"""Util functions."""

import pandas as pd
import cPickle

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


def feature_size(store_name):
    """Number of features."""
    store = pd.io.pytables.HDFStore(store_name)
    a = store.features
    store.close()
    return a.shape[0]


def save_serendipity_dic(y, filename):
    """Save the dictionary with the serendipity values."""
    store = pd.io.pytables.HDFStore(y)
    mat = store.matrix
    store.close()
    n = len(mat.columns)
    ser = 1 - mat.sum(axis=1) / n

    f = open(filename, "w")
    cPickle.dump(ser.to_dict(), f, protocol=2)
    f.close()


def get_serendipity_val(dic, key):
    """Return the serendipity of a RNA."""
    # The key was in the training set
    try:
        return dic[key]
    # The key wasn't in the training set, then the serendipity is 1
    except KeyError:
        return 1.
