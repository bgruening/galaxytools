"""
Z_RandomOverSampler
"""

import imblearn
import numpy as np

from collections import Counter
from imblearn.over_sampling.base import BaseOverSampler
from imblearn.over_sampling import RandomOverSampler
from imblearn.pipeline import Pipeline as imbPipeline
from imblearn.utils import check_target_type
from scipy import sparse
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.preprocessing.data import _handle_zeros_in_scale
from sklearn.utils import check_array, safe_indexing
from sklearn.utils.fixes import nanpercentile
from sklearn.utils.validation import (check_is_fitted, check_X_y,
                                      FLOAT_DTYPES)


class Z_RandomOverSampler(BaseOverSampler):

    def __init__(self, sampling_strategy='auto',
                 return_indices=False,
                 random_state=None,
                 ratio=None,
                 negative_thres=0,
                 positive_thres=-1):
        super(Z_RandomOverSampler, self).__init__(
            sampling_strategy=sampling_strategy, ratio=ratio)
        self.random_state = random_state
        self.return_indices = return_indices
        self.negative_thres = negative_thres
        self.positive_thres = positive_thres

    @staticmethod
    def _check_X_y(X, y):
        y, binarize_y = check_target_type(y, indicate_one_vs_all=True)
        X, y = check_X_y(X, y, accept_sparse=['csr', 'csc'], dtype=None)
        return X, y, binarize_y

    def _fit_resample(self, X, y):
        n_samples = X.shape[0]

        # convert y to z_score
        y_z = (y - y.mean()) / y.std()

        index0 = np.arange(n_samples)
        index_negative = index0[y_z > self.negative_thres]
        index_positive = index0[y_z <= self.positive_thres]
        index_unclassified = [x for x in index0
                              if x not in index_negative
                              and x not in index_positive]

        y_z[index_negative] = 0
        y_z[index_positive] = 1
        y_z[index_unclassified] = -1

        ros = RandomOverSampler(
            sampling_strategy=self.sampling_strategy,
            random_state=self.random_state,
            ratio=self.ratio)
        _, _ = ros.fit_resample(X, y_z)
        sample_indices = ros.sample_indices_

        print("Before sampler: %s. Total after: %s"
              % (Counter(y_z), sample_indices.shape))

        self.sample_indices_ = np.array(sample_indices)

        if self.return_indices:
            return (safe_indexing(X, sample_indices),
                    safe_indexing(y, sample_indices),
                    sample_indices)
        return (safe_indexing(X, sample_indices),
                safe_indexing(y, sample_indices))


def _get_quantiles(X, quantile_range):
    """
    Calculate column percentiles for 2d array

    Parameters
    ----------
    X : array-like, shape [n_samples, n_features]
    """
    quantiles = []
    for feature_idx in range(X.shape[1]):
        if sparse.issparse(X):
            column_nnz_data = X.data[
                X.indptr[feature_idx]: X.indptr[feature_idx + 1]]
            column_data = np.zeros(shape=X.shape[0], dtype=X.dtype)
            column_data[:len(column_nnz_data)] = column_nnz_data
        else:
            column_data = X[:, feature_idx]
        quantiles.append(nanpercentile(column_data, quantile_range))

    quantiles = np.transpose(quantiles)

    return quantiles


class TDMScaler(BaseEstimator, TransformerMixin):
    """
    Scale features using Training Distribution Matching (TDM) algorithm

    References
    ----------
    .. [1] Thompson JA, Tan J and Greene CS (2016) Cross-platform
           normalization of microarray and RNA-seq data for machine
           learning applications. PeerJ 4, e1621.
    """

    def __init__(self, q_lower=25.0, q_upper=75.0, ):
        self.q_lower = q_lower
        self.q_upper = q_upper

    def fit(self, X, y=None):
        """
        Parameters
        ----------
        X : array-like, shape [n_samples, n_features]
        """
        X = check_array(X, copy=True, estimator=self, dtype=FLOAT_DTYPES,
                        force_all_finite=True)

        if not 0 <= self.q_lower <= self.q_upper <= 100:
            raise ValueError("Invalid quantile parameter values: "
                             "q_lower %s, q_upper: %s"
                             % (str(self.q_lower), str(self.q_upper)))

        # TODO sparse data
        quantiles = nanpercentile(X, (self.q_lower, self.q_upper))
        iqr = quantiles[1] - quantiles[0]

        self.q_lower_ = quantiles[0]
        self.q_upper_ = quantiles[1]
        self.iqr_ = _handle_zeros_in_scale(iqr, copy=False)

        self.max_ = np.nanmax(X)
        self.min_ = np.nanmin(X)

        return self

    def transform(self, X):
        """
        Parameters
        ----------
        X : {array-like, sparse matrix}
            The data used to scale along the specified axis.
        """
        check_is_fitted(self, 'iqr_', 'max_')
        X = check_array(X, copy=True, estimator=self, dtype=FLOAT_DTYPES,
                        force_all_finite=True)

        # TODO sparse data
        train_upper_scale = (self.max_ - self.q_upper_) / self.iqr_
        train_lower_scale = (self.q_lower_ - self.min_) / self.iqr_

        test_quantiles = nanpercentile(X, (self.q_lower, self.q_upper))
        test_iqr = _handle_zeros_in_scale(
            test_quantiles[1] - test_quantiles[0], copy=False)

        test_upper_bound = test_quantiles[1] + train_upper_scale * test_iqr
        test_lower_bound = test_quantiles[0] - train_lower_scale * test_iqr

        test_min = np.nanmin(X)
        if test_lower_bound < test_min:
            test_lower_bound = test_min

        X[X > test_upper_bound] = test_upper_bound
        X[X < test_lower_bound] = test_lower_bound

        X = (X - test_lower_bound) / (test_upper_bound - test_lower_bound)\
            * (self.max_ - self.min_) + self.min_

        return X

    def inverse_transform(self, X):
        """
        Scale the data back to the original state
        """
        raise NotImplementedError("Inverse transformation is not implemented!")
