"""
DyRFE
DyRFECV
MyPipeline
MyimbPipeline
check_feature_importances
"""
import numpy as np

from imblearn import under_sampling, over_sampling, combine
from imblearn.pipeline import Pipeline as imbPipeline
from sklearn import (cluster, compose, decomposition, ensemble,
                     feature_extraction, feature_selection,
                     gaussian_process, kernel_approximation,
                     metrics, model_selection, naive_bayes,
                     neighbors, pipeline, preprocessing,
                     svm, linear_model, tree, discriminant_analysis)

from sklearn.base import BaseEstimator
from sklearn.base import MetaEstimatorMixin, clone, is_classifier
from sklearn.feature_selection.rfe import _rfe_single_fit, RFE, RFECV
from sklearn.model_selection import check_cv
from sklearn.metrics.scorer import check_scoring
from sklearn.utils import check_X_y, safe_indexing, safe_sqr
from sklearn.utils._joblib import Parallel, delayed, effective_n_jobs


class DyRFE(RFE):
    """
    Mainly used with DyRFECV

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a ``fit`` method that provides
        information about feature importance either through a ``coef_``
        attribute or through a ``feature_importances_`` attribute.
    n_features_to_select : int or None (default=None)
        The number of features to select. If `None`, half of the features
        are selected.
    step : int, float or list, optional (default=1)
        If greater than or equal to 1, then ``step`` corresponds to the
        (integer) number of features to remove at each iteration.
        If within (0.0, 1.0), then ``step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.
        If list, a series of steps of features to remove at each iteration.
        Iterations stops when steps finish
    verbose : int, (default=0)
        Controls verbosity of output.

    """
    def __init__(self, estimator, n_features_to_select=None, step=1,
                 verbose=0):
        super(DyRFE, self).__init__(estimator, n_features_to_select,
                                    step, verbose)

    def _fit(self, X, y, step_score=None):

        if type(self.step) is not list:
            return super(DyRFE, self)._fit(X, y, step_score)

        # dynamic step
        X, y = check_X_y(X, y, "csc")
        # Initialization
        n_features = X.shape[1]
        if self.n_features_to_select is None:
            n_features_to_select = n_features // 2
        else:
            n_features_to_select = self.n_features_to_select

        step = []
        for s in self.step:
            if 0.0 < s < 1.0:
                step.append(int(max(1, s * n_features)))
            else:
                step.append(int(s))
            if s <= 0:
                raise ValueError("Step must be >0")

        support_ = np.ones(n_features, dtype=np.bool)
        ranking_ = np.ones(n_features, dtype=np.int)

        if step_score:
            self.scores_ = []

        step_i = 0
        # Elimination
        while np.sum(support_) > n_features_to_select and step_i < len(step):

            # if last step is 1, will keep loop
            if step_i == len(step) - 1 and step[step_i] != 0:
                step.append(step[step_i])

            # Remaining features
            features = np.arange(n_features)[support_]

            # Rank the remaining features
            estimator = clone(self.estimator)
            if self.verbose > 0:
                print("Fitting estimator with %d features." % np.sum(support_))

            estimator.fit(X[:, features], y)

            # Get coefs
            if hasattr(estimator, 'coef_'):
                coefs = estimator.coef_
            else:
                coefs = getattr(estimator, 'feature_importances_', None)
            if coefs is None:
                raise RuntimeError('The classifier does not expose '
                                   '"coef_" or "feature_importances_" '
                                   'attributes')

            # Get ranks
            if coefs.ndim > 1:
                ranks = np.argsort(safe_sqr(coefs).sum(axis=0))
            else:
                ranks = np.argsort(safe_sqr(coefs))

            # for sparse case ranks is matrix
            ranks = np.ravel(ranks)

            # Eliminate the worse features
            threshold =\
                min(step[step_i], np.sum(support_) - n_features_to_select)

            # Compute step score on the previous selection iteration
            # because 'estimator' must use features
            # that have not been eliminated yet
            if step_score:
                self.scores_.append(step_score(estimator, features))
            support_[features[ranks][:threshold]] = False
            ranking_[np.logical_not(support_)] += 1

            step_i += 1

        # Set final attributes
        features = np.arange(n_features)[support_]
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(X[:, features], y)

        # Compute step score when only n_features_to_select features left
        if step_score:
            self.scores_.append(step_score(self.estimator_, features))
        self.n_features_ = support_.sum()
        self.support_ = support_
        self.ranking_ = ranking_

        return self


class DyRFECV(RFECV, MetaEstimatorMixin):
    """
    Compared with RFECV, DyRFECV offers flexiable `step` to eleminate
    features, in the format of list, while RFECV supports only fixed number
    of `step`.

    Parameters
    ----------
    estimator : object
        A supervised learning estimator with a ``fit`` method that provides
        information about feature importance either through a ``coef_``
        attribute or through a ``feature_importances_`` attribute.
    step : int or float, optional (default=1)
        If greater than or equal to 1, then ``step`` corresponds to the
        (integer) number of features to remove at each iteration.
        If within (0.0, 1.0), then ``step`` corresponds to the percentage
        (rounded down) of features to remove at each iteration.
        If list, a series of step to remove at each iteration. iteration stopes
        when finishing all steps
        Note that the last iteration may remove fewer than ``step`` features in
        order to reach ``min_features_to_select``.
    min_features_to_select : int, (default=1)
        The minimum number of features to be selected. This number of features
        will always be scored, even if the difference between the original
        feature count and ``min_features_to_select`` isn't divisible by
        ``step``.
    cv : int, cross-validation generator or an iterable, optional
        Determines the cross-validation splitting strategy.
        Possible inputs for cv are:
        - None, to use the default 3-fold cross-validation,
        - integer, to specify the number of folds.
        - :term:`CV splitter`,
        - An iterable yielding (train, test) splits as arrays of indices.
        For integer/None inputs, if ``y`` is binary or multiclass,
        :class:`sklearn.model_selection.StratifiedKFold` is used. If the
        estimator is a classifier or if ``y`` is neither binary nor multiclass,
        :class:`sklearn.model_selection.KFold` is used.
        Refer :ref:`User Guide <cross_validation>` for the various
        cross-validation strategies that can be used here.
        .. versionchanged:: 0.20
            ``cv`` default value of None will change from 3-fold to 5-fold
            in v0.22.
    scoring : string, callable or None, optional, (default=None)
        A string (see model evaluation documentation) or
        a scorer callable object / function with signature
        ``scorer(estimator, X, y)``.
    verbose : int, (default=0)
        Controls verbosity of output.
    n_jobs : int or None, optional (default=None)
        Number of cores to run in parallel while fitting across folds.
        ``None`` means 1 unless in a :obj:`joblib.parallel_backend` context.
        ``-1`` means using all processors. See :term:`Glossary <n_jobs>`
        for more details.
    """
    def __init__(self, estimator, step=1, min_features_to_select=1, cv='warn',
                 scoring=None, verbose=0, n_jobs=None):
        super(DyRFECV, self).__init__(
            estimator, step=step,
            min_features_to_select=min_features_to_select,
            cv=cv, scoring=scoring, verbose=verbose,
            n_jobs=n_jobs)

    def fit(self, X, y, groups=None):
        """Fit the RFE model and automatically tune the number of selected
           features.
        Parameters
        ----------
        X : {array-like, sparse matrix}, shape = [n_samples, n_features]
            Training vector, where `n_samples` is the number of samples and
            `n_features` is the total number of features.
        y : array-like, shape = [n_samples]
            Target values (integers for classification, real numbers for
            regression).
        groups : array-like, shape = [n_samples], optional
            Group labels for the samples used while splitting the dataset into
            train/test set.
        """
        if type(self.step) is not list:
            return super(DyRFECV, self).fit(X, y, groups)

        X, y = check_X_y(X, y, "csr")

        # Initialization
        cv = check_cv(self.cv, y, is_classifier(self.estimator))
        scorer = check_scoring(self.estimator, scoring=self.scoring)
        n_features = X.shape[1]

        step = []
        for s in self.step:
            if 0.0 < s < 1.0:
                step.append(int(max(1, s * n_features)))
            else:
                step.append(int(s))
            if s <= 0:
                raise ValueError("Step must be >0")

        # Build an RFE object, which will evaluate and score each possible
        # feature count, down to self.min_features_to_select
        rfe = DyRFE(estimator=self.estimator,
                    n_features_to_select=self.min_features_to_select,
                    step=self.step, verbose=self.verbose)

        # Determine the number of subsets of features by fitting across
        # the train folds and choosing the "features_to_select" parameter
        # that gives the least averaged error across all folds.

        # Note that joblib raises a non-picklable error for bound methods
        # even if n_jobs is set to 1 with the default multiprocessing
        # backend.
        # This branching is done so that to
        # make sure that user code that sets n_jobs to 1
        # and provides bound methods as scorers is not broken with the
        # addition of n_jobs parameter in version 0.18.

        if effective_n_jobs(self.n_jobs) == 1:
            parallel, func = list, _rfe_single_fit
        else:
            parallel = Parallel(n_jobs=self.n_jobs)
            func = delayed(_rfe_single_fit)

        scores = parallel(
            func(rfe, self.estimator, X, y, train, test, scorer)
            for train, test in cv.split(X, y, groups))

        scores = np.sum(scores, axis=0)
        diff = int(scores.shape[0]) - len(step)
        if diff > 0:
            step = np.r_[step, [step[-1]] * diff]
        scores_rev = scores[::-1]
        argmax_idx = len(scores) - np.argmax(scores_rev) - 1
        n_features_to_select = max(
            n_features - sum(step[:argmax_idx]),
            self.min_features_to_select)

        # Re-execute an elimination with best_k over the whole set
        rfe = DyRFE(estimator=self.estimator,
                    n_features_to_select=n_features_to_select, step=self.step,
                    verbose=self.verbose)

        rfe.fit(X, y)

        # Set final attributes
        self.support_ = rfe.support_
        self.n_features_ = rfe.n_features_
        self.ranking_ = rfe.ranking_
        self.estimator_ = clone(self.estimator)
        self.estimator_.fit(self.transform(X), y)

        # Fixing a normalization error, n is equal to get_n_splits(X, y) - 1
        # here, the scores are normalized by get_n_splits(X, y)
        self.grid_scores_ = scores[::-1] / cv.get_n_splits(X, y, groups)
        return self


class MyPipeline(pipeline.Pipeline):
    """
    Extend pipeline object to have feature_importances_ attribute
    """
    def fit(self, X, y=None, **fit_params):
        super(MyPipeline, self).fit(X, y, **fit_params)
        estimator = self.steps[-1][-1]
        if hasattr(estimator, 'coef_'):
            coefs = estimator.coef_
        else:
            coefs = getattr(estimator, 'feature_importances_', None)
        if coefs is None:
            raise RuntimeError('The estimator in the pipeline does not expose '
                               '"coef_" or "feature_importances_" '
                               'attributes')
        self.feature_importances_ = coefs
        return self


class MyimbPipeline(imbPipeline):
    """
    Extend imblance pipeline object to have feature_importances_ attribute
    """
    def fit(self, X, y=None, **fit_params):
        super(MyimbPipeline, self).fit(X, y, **fit_params)
        estimator = self.steps[-1][-1]
        if hasattr(estimator, 'coef_'):
            coefs = estimator.coef_
        else:
            coefs = getattr(estimator, 'feature_importances_', None)
        if coefs is None:
            raise RuntimeError('The estimator in the pipeline does not expose '
                               '"coef_" or "feature_importances_" '
                               'attributes')
        self.feature_importances_ = coefs
        return self


def check_feature_importances(estimator):
    """
    For pipeline object which has no feature_importances_ property,
    this function returns the same comfigured pipeline object with
    attached the last estimator's feature_importances_.
    """
    if estimator.__class__.__module__ == 'sklearn.pipeline':
        pipeline_steps = estimator.get_params()['steps']
        estimator = MyPipeline(pipeline_steps)
    elif estimator.__class__.__module__ == 'imblearn.pipeline':
        pipeline_steps = estimator.get_params()['steps']
        estimator = MyimbPipeline(pipeline_steps)
    else:
        return estimator
