#!/usr/bin/env python
"""Compute scipy.stats tests and descriptive statistics for tabular data.

The selected columns of every row of the input file are treated as one
sample (or several samples, separated by ``;`` for the multi-sample
tests) and the result is appended to the row, or, in per-column mode,
each selected column (all of its values) forms one sample. The computed
result is written to the output file, optionally with a header row.
"""

import argparse

import numpy as np
from scipy import stats

# Must match @TOOL_VERSION@ in macros.xml
__version__ = "1.16.3"

RESULT_LABELS = {
    "anderson": ["statistic", "critical_values", "significance_levels"],
    "ansari": ["statistic", "pvalue"],
    "bartlett": ["statistic", "pvalue"],
    "bayes_mvs": ["mean", "variance", "std_dev"],
    "binom_test": ["pvalue"],
    "binned_statistic": ["statistic", "bin_edges", "binnumber"],
    "bootstrap": ["ci_low", "ci_high", "standard_error"],
    "boxcox": ["transformed"],
    "boxcox_llf": ["llf"],
    "boxcox_normmax": ["maxlog"],
    "chisquare": ["statistic", "pvalue"],
    "chi2_contingency": ["chi2", "pvalue", "dof", "expected"],
    "combine_pvalues": ["statistic", "pvalue"],
    "cumfreq": ["cumfreq", "lowerlimit", "binsize", "extrapoints"],
    "describe": ["size", "min", "max", "mean", "variance", "skewness", "kurtosis"],
    "entropy": ["entropy"],
    "f_oneway": ["statistic", "pvalue"],
    "fligner": ["statistic", "pvalue"],
    "friedmanchisquare": ["statistic", "pvalue"],
    "gmean": ["gmean"],
    "histogram": ["counts", "lowerlimit", "binsize", "extrapoints"],
    "histogram2": ["counts"],
    "hmean": ["hmean"],
    "itemfreq": ["values", "counts"],
    "kendalltau": ["correlation", "pvalue"],
    "kruskal": ["statistic", "pvalue"],
    "kstest": ["statistic", "pvalue"],
    "ks_2samp": ["statistic", "pvalue"],
    "kurtosis": ["kurtosis"],
    "kurtosistest": ["statistic", "pvalue"],
    "levene": ["statistic", "pvalue"],
    "linregress": ["slope", "intercept", "rvalue", "pvalue", "stderr"],
    "mannwhitneyu": ["statistic", "pvalue"],
    "median_test": ["statistic", "pvalue", "median", "table"],
    "mode": ["mode", "count"],
    "moment": ["moment"],
    "mood": ["statistic", "pvalue"],
    "nanmean": ["nanmean"],
    "nanmedian": ["nanmedian"],
    "nanstd": ["nanstd"],
    "normaltest": ["statistic", "pvalue"],
    "obrientransform": ["transformed"],
    "pearsonr": ["correlation", "pvalue"],
    "percentileofscore": ["percentile"],
    "permutation_test": ["statistic", "pvalue"],
    "pointbiserialr": ["correlation", "pvalue"],
    "power_divergence": ["statistic", "pvalue"],
    "rankdata": ["ranks"],
    "ranksums": ["statistic", "pvalue"],
    "relfreq": ["frequency", "lowerlimit", "binsize", "extrapoints"],
    "scoreatpercentile": ["score"],
    "sem": ["sem"],
    "shapiro": ["statistic", "pvalue"],
    "signaltonoise": ["signaltonoise"],
    "sigmaclip": ["clipped", "lower", "upper"],
    "skew": ["skewness"],
    "skewtest": ["statistic", "pvalue"],
    "spearmanr": ["correlation", "pvalue"],
    "theilslopes": ["medslope", "medintercept", "low_slope", "high_slope"],
    "threshold": ["thresholded"],
    "tiecorrect": ["tie_correction"],
    "trim1": ["trimmed"],
    "trimboth": ["trimmed"],
    "tmean": ["tmean"],
    "tmax": ["tmax"],
    "tmin": ["tmin"],
    "tstd": ["tstd"],
    "tsem": ["tsem"],
    "tvar": ["tvar"],
    "ttest_1samp": ["t_statistic", "pvalue"],
    "ttest_ind": ["t_statistic", "pvalue"],
    "ttest_rel": ["t_statistic", "pvalue"],
    "variation": ["variation"],
    "wilcoxon": ["statistic", "pvalue"],
    "zmap": ["zscore"],
    "zscore": ["zscore"],
}

BOOTSTRAP_STATISTICS = {
    "mean": np.mean,
    "median": np.median,
    "std": np.std,
    "var": np.var,
    "min": np.min,
    "max": np.max,
    "sum": np.sum,
    "count": np.size,
}

PERMUTATION_STATISTICS = {
    "mean": lambda x, y, axis=-1: np.mean(x, axis=axis) - np.mean(y, axis=axis),
    "median": lambda x, y, axis=-1: np.median(x, axis=axis) - np.median(y, axis=axis),
}


def parse_columns(spec):
    if spec is None:
        return []
    return [int(col) for col in str(spec).split(",") if col.strip() != ""]


def parse_float(cell):
    cell = cell.strip()
    if not cell:
        return float("nan")
    return float(cell)


def stringify(value):
    if isinstance(value, np.ndarray):
        return ",".join(stringify(v) for v in value.tolist())
    if isinstance(value, np.generic):
        return str(value.item())
    if isinstance(value, (tuple, list)):
        return "(" + ",".join(stringify(v) for v in value) + ")"
    return str(value)


def signaltonoise(a, ddof=0):
    a = np.asarray(a, dtype=float)
    mean = a.mean()
    sd = a.std(ddof=ddof)
    if sd == 0:
        return 0.0
    return abs(mean) / sd


def threshold(a, threshmin=None, threshmax=None, newval=0.0):
    """Replacement for the scipy.stats.threshold function removed in
    scipy 1.12: values below threshmin or above threshmax are replaced
    by newval."""
    a = np.array(a, dtype=float)
    mask = np.ones(len(a), dtype=bool)
    if threshmin is not None:
        mask &= a >= threshmin
    if threshmax is not None:
        mask &= a <= threshmax
    a[~mask] = newval
    return a


def histogram(a, numbins=10, defaultlimits=None):
    """Replacement for the scipy.stats.histogram function removed in
    scipy 1.12."""
    a = np.asarray(a, dtype=float)
    if defaultlimits is None:
        low, high = a.min(), a.max()
    else:
        low, high = defaultlimits
    counts, _edges = np.histogram(a, bins=numbins, range=(low, high))
    binsize = (high - low) / numbins
    extrapoints = len(a) - int(counts.sum())
    return counts.astype(float), low, binsize, extrapoints


def histogram2(a, bins):
    """Replacement for the scipy.stats.histogram2 function removed in
    scipy 1.12."""
    n = np.searchsorted(np.sort(a), bins)
    n = np.concatenate([n, [len(a)]])
    return n[1:] - n[:-1]


def apply_nan_policy(values, policy):
    values = np.asarray(values, dtype=float)
    if policy == "omit":
        return [v for v in values.tolist() if not np.isnan(v)]
    if policy == "raise" and np.isnan(values).any():
        raise ValueError("sample contains NaN values")
    return values.tolist()


def bonferroni(pvalues):
    p = np.asarray(pvalues, dtype=float)
    return np.minimum(p * len(p), 1.0)


def benjamini_hochberg(pvalues):
    p = np.asarray(pvalues, dtype=float)
    n = len(p)
    order = np.argsort(p, kind="stable")
    ranked = p[order] * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(adjusted, 0.0, 1.0)
    return out


def make_rng(seed):
    return np.random.default_rng(seed) if seed is not None else None


def compute(test_id, args, sample_one, sample_two, multi_samples):
    """Run the selected test on one sample and return the result fields."""

    def floats(values):
        return list(map(float, values))

    if test_id == "describe":
        size, min_max, mean, uv, bs, bk = stats.describe(floats(sample_one))
        return [size, min_max, mean, uv, bs, bk]
    if test_id == "mode":
        vals, counts = stats.mode(floats(sample_one))
        return [vals, counts]
    if test_id == "nanmean":
        return [np.nanmean(floats(sample_one))]
    if test_id == "nanmedian":
        return [np.nanmedian(floats(sample_one))]
    if test_id == "nanstd":
        ddof = 0 if args.bias else 1
        return [np.nanstd(floats(sample_one), ddof=ddof)]
    if test_id == "kurtosistest":
        z_value, p_value = stats.kurtosistest(floats(sample_one))
        return [z_value, p_value]
    if test_id == "variation":
        return [stats.variation(floats(sample_one))]
    if test_id == "itemfreq":
        freq = np.unique(floats(sample_one), return_counts=True)
        return [",".join(map(str, i)) for i in freq]
    if test_id == "boxcox_llf":
        return [stats.boxcox_llf(args.imbda, floats(sample_one))]
    if test_id == "tiecorrect":
        return [stats.tiecorrect(floats(sample_one))]
    if test_id == "rankdata":
        return [stats.rankdata(floats(sample_one), method=args.md)]
    if test_id == "anderson":
        A2, critical, sig = stats.anderson(floats(sample_one), dist=args.dist)
        return [A2, critical, sig]
    if test_id == "binom_test":
        x = [int(round(v)) for v in sample_one]
        if len(x) == 2:
            successes, n = x[0], x[0] + x[1]
        elif len(x) == 1:
            if args.n is None:
                raise ValueError(
                    "For a single binom_test column the number of trials must be set"
                )
            successes, n = x[0], args.n
        else:
            raise ValueError("binom_test needs one or two columns")
        return [stats.binomtest(successes, n, p=args.p).pvalue]
    if test_id == "gmean":
        return [stats.gmean(floats(sample_one), dtype=args.dtype)]
    if test_id == "hmean":
        return [stats.hmean(floats(sample_one), dtype=args.dtype)]
    if test_id == "kurtosis":
        return [
            stats.kurtosis(
                floats(sample_one),
                axis=args.axis,
                fisher=args.fisher,
                bias=args.bias,
            )
        ]
    if test_id == "moment":
        order = args.n if args.n is not None else 1
        return [stats.moment(floats(sample_one), order=order)]
    if test_id == "normaltest":
        k2, p_value = stats.normaltest(floats(sample_one))
        return [k2, p_value]
    if test_id == "skew":
        return [stats.skew(floats(sample_one), bias=args.bias)]
    if test_id == "skewtest":
        z_value, p_value = stats.skewtest(floats(sample_one))
        return [z_value, p_value]
    if test_id == "sem":
        return [stats.sem(floats(sample_one), ddof=args.ddof)]
    if test_id == "zscore":
        return list(stats.zscore(floats(sample_one), ddof=args.ddof))
    if test_id == "signaltonoise":
        return [signaltonoise(floats(sample_one), ddof=args.ddof)]
    if test_id == "percentileofscore":
        return [
            stats.percentileofscore(
                floats(sample_one), score=args.score, kind=args.kind
            )
        ]
    if test_id == "bayes_mvs":
        c_mean, c_var, c_std = stats.bayes_mvs(floats(sample_one), alpha=args.alpha)
        return [c_mean, c_var, c_std]
    if test_id == "sigmaclip":
        low = args.m if args.m is not None else 4.0
        high = args.n if args.n is not None else 4.0
        c, c_low, c_up = stats.sigmaclip(floats(sample_one), low=low, high=high)
        return [c, c_low, c_up]
    if test_id == "kstest":
        d, p_value = stats.kstest(
            floats(sample_one),
            cdf=args.cdf,
            N=args.ni,
            alternative=args.alternative,
            method=args.mode if args.mode else "auto",
        )
        return [d, p_value]
    if test_id == "chi2_contingency":
        chi2, p, dof, ex = stats.chi2_contingency(
            floats(sample_one),
            correction=args.correction,
            lambda_=args.lambda_,
        )
        return [chi2, p, dof, ex]
    if test_id == "tmean":
        if args.mf is None and args.nf is None:
            return [stats.tmean(floats(sample_one))]
        lower = args.mf if args.mf is not None else -np.inf
        upper = args.nf if args.nf is not None else np.inf
        return [
            stats.tmean(
                floats(sample_one),
                (lower, upper),
                (args.inclusive1, args.inclusive2),
            )
        ]
    if test_id == "tmin":
        return [
            stats.tmin(
                floats(sample_one),
                lowerlimit=args.mf,
                inclusive=args.inclusive,
            )
        ]
    if test_id == "tmax":
        return [
            stats.tmax(
                floats(sample_one),
                upperlimit=args.nf,
                inclusive=args.inclusive,
            )
        ]
    if test_id == "tvar":
        if args.mf is None and args.nf is None:
            return [stats.tvar(floats(sample_one))]
        lower = args.mf if args.mf is not None else -np.inf
        upper = args.nf if args.nf is not None else np.inf
        return [
            stats.tvar(
                floats(sample_one),
                (lower, upper),
                (args.inclusive1, args.inclusive2),
            )
        ]
    if test_id == "tstd":
        if args.mf is None and args.nf is None:
            return [stats.tstd(floats(sample_one))]
        lower = args.mf if args.mf is not None else -np.inf
        upper = args.nf if args.nf is not None else np.inf
        return [
            stats.tstd(
                floats(sample_one),
                (lower, upper),
                (args.inclusive1, args.inclusive2),
            )
        ]
    if test_id == "tsem":
        if args.mf is None and args.nf is None:
            return [stats.tsem(floats(sample_one))]
        lower = args.mf if args.mf is not None else -np.inf
        upper = args.nf if args.nf is not None else np.inf
        return [
            stats.tsem(
                floats(sample_one),
                (lower, upper),
                (args.inclusive1, args.inclusive2),
            )
        ]
    if test_id == "scoreatpercentile":
        limit = (
            (args.mf, args.nf)
            if args.mf is not None and args.nf is not None
            else ()
        )
        scores = stats.scoreatpercentile(
            floats(sample_one),
            floats(sample_two),
            limit,
            interpolation_method=args.interpolation,
        )
        return list(scores)
    if test_id == "relfreq":
        limits = (
            (args.mf, args.nf)
            if args.mf is not None and args.nf is not None
            else None
        )
        rel, low_range, binsize, ex = stats.relfreq(
            floats(sample_one), args.b, limits
        )
        return list(rel) + [low_range, binsize, ex]
    if test_id == "binned_statistic":
        limits = (
            (args.mf, args.nf)
            if args.mf is not None and args.nf is not None
            else None
        )
        st, b_edge, b_n = stats.binned_statistic(
            floats(sample_one),
            floats(sample_two),
            statistic=args.statistic,
            bins=args.b,
            range=limits,
        )
        return [st, b_edge, b_n]
    if test_id == "threshold":
        return list(
            threshold(
                floats(sample_one),
                threshmin=args.mf,
                threshmax=args.nf,
                newval=args.new,
            )
        )
    if test_id == "trimboth":
        return list(
            stats.trimboth(
                floats(sample_one), proportiontocut=args.proportiontocut
            )
        )
    if test_id == "trim1":
        return list(
            stats.trim1(
                floats(sample_one),
                proportiontocut=args.proportiontocut,
                tail=args.tail,
            )
        )
    if test_id == "histogram":
        limits = (
            (args.mf, args.nf)
            if args.mf is not None and args.nf is not None
            else None
        )
        hi, low_range, binsize, ex = histogram(floats(sample_one), args.b, limits)
        return [hi, low_range, binsize, ex]
    if test_id == "cumfreq":
        limits = (
            (args.mf, args.nf)
            if args.mf is not None and args.nf is not None
            else None
        )
        cum, low_range, binsize, ex = stats.cumfreq(
            floats(sample_one), args.b, limits
        )
        return [cum, low_range, binsize, ex]
    if test_id == "boxcox_normmax":
        if args.mf is not None and args.nf is not None:
            ma = stats.boxcox_normmax(
                floats(sample_one), (args.mf, args.nf), method=args.method
            )
        else:
            ma = stats.boxcox_normmax(floats(sample_one), method=args.method)
        return [ma]
    if test_id == "boxcox":
        if args.imbda:
            return [stats.boxcox(floats(sample_one), lmbda=args.imbda)]
        box, ma, ci = stats.boxcox(floats(sample_one), alpha=args.alpha)
        return [box, ma, ci]
    if test_id == "histogram2":
        return list(histogram2(floats(sample_one), floats(sample_two)))
    if test_id == "ranksums":
        z_statistic, p_value = stats.ranksums(
            floats(sample_one), floats(sample_two)
        )
        return [z_statistic, p_value]
    if test_id == "ttest_1samp":
        if len(sample_two) != 1:
            raise ValueError(
                "ttest_1samp needs exactly one column as popmean"
            )
        result = stats.ttest_1samp(floats(sample_one), sample_two[0])
        results = [result.statistic, result.pvalue]
        if args.confidence_level is not None:
            ci = result.confidence_interval(confidence_level=args.confidence_level)
            results += [ci.low, ci.high]
        return results
    if test_id == "ansari":
        AB, p_value = stats.ansari(floats(sample_one), floats(sample_two))
        return [AB, p_value]
    if test_id == "linregress":
        slope, intercept, r_value, p_value, stderr = stats.linregress(
            floats(sample_one), floats(sample_two)
        )
        return [slope, intercept, r_value, p_value, stderr]
    if test_id == "pearsonr":
        cor, p_value = stats.pearsonr(floats(sample_one), floats(sample_two))
        return [cor, p_value]
    if test_id == "pointbiserialr":
        r, p_value = stats.pointbiserialr(
            floats(sample_one), floats(sample_two)
        )
        return [r, p_value]
    if test_id == "ks_2samp":
        d, p_value = stats.ks_2samp(floats(sample_one), floats(sample_two))
        return [d, p_value]
    if test_id == "mannwhitneyu":
        mw_stats_u, p_value = stats.mannwhitneyu(
            floats(sample_one),
            floats(sample_two),
            use_continuity=args.mwu_use_continuity,
        )
        return [mw_stats_u, p_value]
    if test_id == "zmap":
        return list(
            stats.zmap(floats(sample_one), floats(sample_two), ddof=args.ddof)
        )
    if test_id == "ttest_ind":
        mw_stats_u, p_value = stats.ttest_ind(
            floats(sample_one),
            floats(sample_two),
            equal_var=args.equal_var,
        )
        return [mw_stats_u, p_value]
    if test_id == "ttest_rel":
        t, prob = stats.ttest_rel(
            floats(sample_one), floats(sample_two), axis=args.axis
        )
        return [t, prob]
    if test_id == "mood":
        z, p_value = stats.mood(
            floats(sample_one), floats(sample_two), axis=args.axis
        )
        return [z, p_value]
    if test_id == "shapiro":
        W, p_value = stats.shapiro(floats(sample_one))
        return [W, p_value]
    if test_id == "kendalltau":
        k, p_value = stats.kendalltau(floats(sample_one), floats(sample_two))
        return [k, p_value]
    if test_id == "entropy":
        qk = floats(sample_two) if sample_two else None
        return [stats.entropy(floats(sample_one), qk, base=args.base)]
    if test_id == "spearmanr":
        if not sample_two:
            raise ValueError("spearmanr needs two samples")
        rho, p_value = stats.spearmanr(
            floats(sample_one), floats(sample_two)
        )
        return [rho, p_value]
    if test_id == "wilcoxon":
        if sample_two:
            T, p_value = stats.wilcoxon(
                floats(sample_one),
                floats(sample_two),
                zero_method=args.zero_method,
                correction=args.correction,
            )
        else:
            T, p_value = stats.wilcoxon(
                floats(sample_one),
                zero_method=args.zero_method,
                correction=args.correction,
            )
        return [T, p_value]
    if test_id == "chisquare":
        if sample_two:
            rho, p_value = stats.chisquare(
                floats(sample_one), floats(sample_two), ddof=args.ddof
            )
        else:
            rho, p_value = stats.chisquare(
                floats(sample_one), ddof=args.ddof
            )
        return [rho, p_value]
    if test_id == "power_divergence":
        if sample_two:
            stat, p_value = stats.power_divergence(
                floats(sample_one),
                floats(sample_two),
                ddof=args.ddof,
                lambda_=args.lambda_,
            )
        else:
            stat, p_value = stats.power_divergence(
                floats(sample_one), ddof=args.ddof, lambda_=args.lambda_
            )
        return [stat, p_value]
    if test_id == "theilslopes":
        if sample_two:
            mpe, met, lo, up = stats.theilslopes(
                floats(sample_one), floats(sample_two), alpha=args.alpha
            )
        else:
            mpe, met, lo, up = stats.theilslopes(
                floats(sample_one), alpha=args.alpha
            )
        return [mpe, met, lo, up]
    if test_id == "combine_pvalues":
        if sample_two:
            stat, p_value = stats.combine_pvalues(
                floats(sample_one),
                method=args.med,
                weights=floats(sample_two),
            )
        else:
            stat, p_value = stats.combine_pvalues(
                floats(sample_one), method=args.med
            )
        return [stat, p_value]
    if test_id == "bootstrap":
        statistic = BOOTSTRAP_STATISTICS.get(args.statistic)
        if statistic is None:
            raise ValueError(
                "unsupported statistic for bootstrap: %s" % args.statistic
            )
        res = stats.bootstrap(
            (np.asarray(floats(sample_one)),),
            statistic,
            n_resamples=args.n_resamples,
            confidence_level=args.confidence_level,
            rng=make_rng(args.seed),
        )
        return [res.confidence_interval.low, res.confidence_interval.high, res.standard_error]
    if test_id == "permutation_test":
        statistic = PERMUTATION_STATISTICS.get(args.statistic)
        if statistic is None:
            raise ValueError(
                "unsupported statistic for permutation_test: %s" % args.statistic
            )
        res = stats.permutation_test(
            (np.asarray(floats(sample_one)), np.asarray(floats(sample_two))),
            statistic,
            permutation_type=args.permutation_type,
            n_resamples=args.n_resamples,
            alternative=args.alternative,
            rng=make_rng(args.seed),
        )
        return [res.statistic, res.pvalue]
    if test_id == "obrientransform":
        ob = stats.obrientransform(*multi_samples)
        return [",".join(map(str, i)) for i in ob]
    if test_id == "f_oneway":
        f_value, p_value = stats.f_oneway(*multi_samples)
        return [f_value, p_value]
    if test_id == "kruskal":
        h, p_value = stats.kruskal(*multi_samples)
        return [h, p_value]
    if test_id == "friedmanchisquare":
        fr, p_value = stats.friedmanchisquare(*multi_samples)
        return [fr, p_value]
    if test_id == "fligner":
        xsq, p_value = stats.fligner(
            *multi_samples,
            center=args.center,
            proportiontocut=args.proportiontocut,
        )
        return [xsq, p_value]
    if test_id == "bartlett":
        T, p_value = stats.bartlett(*multi_samples)
        return [T, p_value]
    if test_id == "levene":
        w, p_value = stats.levene(
            *multi_samples,
            center=args.center,
            proportiontocut=args.proportiontocut,
        )
        return [w, p_value]
    if test_id == "median_test":
        stat, p_value, m, table = stats.median_test(
            *multi_samples,
            ties=args.ties,
            correction=args.correction,
            lambda_=args.lambda_,
        )
        return [stat, p_value, m] + [",".join(map(str, i)) for i in table]
    raise ValueError("Unknown test_id: %s" % test_id)


def build_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True, help="Tabular file.")
    parser.add_argument(
        "-o", "--outfile", required=True, help="Path to the output file."
    )
    parser.add_argument("--test_id", help="statistical test method")
    parser.add_argument(
        "--test_scope",
        choices=("per_row", "per_column"),
        default="per_row",
        help="Whether the selected columns of every row form one sample "
        "(per_row) or every selected column over all rows forms one sample "
        "(per_column)",
    )
    parser.add_argument(
        "--has_header",
        action="store_true",
        default=False,
        help="The first line of the input file is a header and is skipped",
    )
    parser.add_argument(
        "--nan_policy",
        choices=("propagate", "omit", "raise"),
        default="propagate",
        help="How to handle NaN values in the samples",
    )
    parser.add_argument(
        "--output_mode",
        choices=("append", "results_only"),
        default="append",
        help="Whether to append the results to the input columns or to "
        "write only the results",
    )
    parser.add_argument(
        "--include_header",
        action="store_true",
        default=False,
        help="Write a header row describing the result columns",
    )
    parser.add_argument(
        "--pvalue_correction",
        choices=("none", "bonferroni", "fdr"),
        default="none",
        help="Multiple testing correction applied to the p-value column",
    )
    parser.add_argument("--sample_one_cols", help="Columns of sample one")
    parser.add_argument("--sample_two_cols", help="Columns of sample two")
    parser.add_argument(
        "--sample_cols", help="Columns of several samples, separated arrays using ;"
    )
    parser.add_argument(
        "--n_resamples",
        type=int,
        default=9999,
        help="Number of resamples used by bootstrap and permutation_test",
    )
    parser.add_argument(
        "--confidence_level",
        type=float,
        default=None,
        help="Confidence level of the confidence interval",
    )
    parser.add_argument(
        "--permutation_type",
        choices=("independent", "samples", "pairwise"),
        default="independent",
        help="Permutation type used by permutation_test",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=None,
        help="Random seed used by bootstrap and permutation_test",
    )
    parser.add_argument(
        "--mwu_use_continuity",
        action="store_true",
        default=False,
        help="Whether a continuity correction (1/2.) should be taken into account.",
    )
    parser.add_argument(
        "--equal_var",
        action="store_true",
        default=False,
        help="If set perform a standard independent 2 sample test that assumes "
        "equal population variances. If not set, perform Welch's t-test, which "
        "does not assume equal population variance.",
    )
    parser.add_argument(
        "--fisher",
        action="store_true",
        default=False,
        help="if true then Fisher definition is used",
    )
    parser.add_argument(
        "--bias",
        action="store_true",
        default=False,
        help="if false,then the calculations are corrected for statistical bias",
    )
    parser.add_argument(
        "--inclusive1",
        action="store_true",
        default=False,
        help="if false,lower_limit will be ignored",
    )
    parser.add_argument(
        "--inclusive2",
        action="store_true",
        default=False,
        help="if false,higher_limit will be ignored",
    )
    parser.add_argument(
        "--inclusive",
        action="store_true",
        default=False,
        help="if false,limit will be ignored",
    )
    parser.add_argument(
        "--printextras",
        action="store_true",
        default=False,
        help="If True, if there are extra points a warning is raised saying "
        "how many of those points there are",
    )
    parser.add_argument(
        "--correction",
        action="store_true",
        default=False,
        help="continuity correction ",
    )
    parser.add_argument(
        "--axis",
        type=int,
        default=0,
        help="Axis along which to operate, 0 for one-dimensional data",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=None,
        help="Number used by the selected test: order of the moment (moment), "
        "number of trials (binom_test) or upper bound factor (sigmaclip)",
    )
    parser.add_argument(
        "--b", type=int, default=10, help="The number of bins to use for the histogram"
    )
    parser.add_argument(
        "--ni", type=int, default=20, help="Sample size used by the kstest"
    )
    parser.add_argument(
        "--ddof", type=int, default=0, help="Degrees of freedom correction"
    )
    parser.add_argument(
        "--score",
        type=float,
        default=0.0,
        help="Score that is compared to the elements in a.",
    )
    parser.add_argument("--m", type=float, default=None, help="limits")
    parser.add_argument("--mf", type=float, default=None, help="lower limit")
    parser.add_argument("--nf", type=float, default=None, help="higher limit")
    parser.add_argument(
        "--p",
        type=float,
        default=0.5,
        help="The hypothesized probability of success. 0 <= p <= 1. "
        "The default value is p = 0.5",
    )
    parser.add_argument("--alpha", type=float, default=0.9, help="probability")
    parser.add_argument(
        "--new",
        type=float,
        default=0.0,
        help="Value to put in place of values in a outside of bounds",
    )
    parser.add_argument(
        "--proportiontocut",
        type=float,
        default=0.0,
        help="Proportion (in range 0-1) of total data set to trim of each end.",
    )
    parser.add_argument(
        "--lambda_",
        type=float,
        default=1.0,
        help="lambda_ gives the power in the Cressie-Read power divergence statistic",
    )
    parser.add_argument(
        "--imbda",
        type=float,
        default=0.0,
        help="If lmbda is not None, do the transformation for that value. "
        "If lmbda is None, find the lambda that maximizes the log-likelihood "
        "function and return it as the second output argument.",
    )
    parser.add_argument(
        "--base",
        type=float,
        default=None,
        help="The logarithmic base to use, defaults to e",
    )
    parser.add_argument("--dtype", default=None, help="dtype")
    parser.add_argument("--med", help="med")
    parser.add_argument("--cdf", default="norm", help="cdf")
    parser.add_argument("--zero_method", help="zero_method options")
    parser.add_argument("--dist", help="dist options")
    parser.add_argument("--ties", help="ties options")
    parser.add_argument("--alternative", help="alternative options")
    parser.add_argument("--mode", help="mode options")
    parser.add_argument("--method", help="method options")
    parser.add_argument("--md", help="md options")
    parser.add_argument("--center", help="center options")
    parser.add_argument("--kind", help="kind options")
    parser.add_argument("--tail", help="tail options")
    parser.add_argument("--interpolation", help="interpolation options")
    parser.add_argument("--statistic", help="statistic options")
    return parser


def column_name(header_fields, index):
    if header_fields is not None and index <= len(header_fields):
        return header_fields[index - 1]
    return "column_%d" % index


def run_per_column(test_id, args, data_lines, header_fields):
    """Compute the test for every selected column over all of its values."""
    sample_one_cols = parse_columns(args.sample_one_cols)
    sample_two_cols = parse_columns(args.sample_two_cols)
    multi_samples = [
        parse_columns(group) for group in (args.sample_cols or "").split(";") if group
    ]
    rows = []
    if multi_samples:
        groups = []
        for group in multi_samples:
            values = [
                parse_float(line.split("\t")[c - 1])
                for line in data_lines
                for c in group
            ]
            groups.append(apply_nan_policy(values, args.nan_policy))
        name = ";".join(
            ",".join(column_name(header_fields, c) for c in group)
            for group in multi_samples
        )
        try:
            results = compute(test_id, args, groups[0], groups[1] if len(groups) > 1 else [], groups)
        except Exception:
            raise ValueError("failed to compute %s" % test_id)
        rows.append((name, results))
    else:
        if len(sample_two_cols) == len(sample_one_cols):
            pairs = list(zip(sample_one_cols, sample_two_cols))
        elif len(sample_two_cols) == 1:
            pairs = [(c, sample_two_cols[0]) for c in sample_one_cols]
        else:
            pairs = [(c, None) for c in sample_one_cols]
        for one_col, two_col in pairs:
            one_values = apply_nan_policy(
                [parse_float(line.split("\t")[one_col - 1]) for line in data_lines],
                args.nan_policy,
            )
            two_values = []
            if two_col is not None:
                two_values = apply_nan_policy(
                    [parse_float(line.split("\t")[two_col - 1]) for line in data_lines],
                    args.nan_policy,
                )
            try:
                results = compute(test_id, args, one_values, two_values, [])
            except Exception:
                raise ValueError(
                    "failed to compute %s for column %d" % (test_id, one_col)
                )
            rows.append((column_name(header_fields, one_col), results))
    return rows


def main():
    args = build_parser().parse_args()
    test_id = (args.test_id or "").strip()
    if test_id not in RESULT_LABELS:
        raise ValueError("Unknown test_id: %s" % test_id)
    labels = list(RESULT_LABELS[test_id])
    if test_id == "ttest_1samp" and args.confidence_level is not None:
        labels += ["ci_low", "ci_high"]
    correction = args.pvalue_correction
    pvalue_index = (
        labels.index("pvalue")
        if correction != "none" and "pvalue" in labels
        else None
    )

    with open(args.infile) as fin:
        all_lines = [line for line in fin.read().splitlines() if line.strip()]
    header_fields = None
    if args.has_header:
        if not all_lines:
            raise ValueError("the input file is empty")
        header_fields = all_lines[0].split("\t")
        data_lines = all_lines[1:]
    else:
        data_lines = all_lines
    if not data_lines:
        raise ValueError("the input file contains no data rows")

    with open(args.outfile, "w") as fout:

        def result_fields(results):
            return [stringify(result) for result in results]

        def header_row(fields):
            if args.include_header:
                fout.write("\t".join(fields) + "\n")

        if args.test_scope == "per_column":
            rows = run_per_column(test_id, args, data_lines, header_fields)
            fields = ["sample"] + labels
            if pvalue_index is not None:
                fields.append("pvalue_adjusted")
            header_row(fields)
            if pvalue_index is None:
                for name, results in rows:
                    fout.write("\t".join([name] + result_fields(results)) + "\n")
            else:
                pvalues = [results[pvalue_index] for _, results in rows]
                if correction == "bonferroni":
                    adjusted = bonferroni(pvalues)
                else:
                    adjusted = benjamini_hochberg(pvalues)
                for (name, results), adj in zip(rows, adjusted):
                    fout.write(
                        "\t".join([name] + result_fields(results) + [stringify(adj)])
                        + "\n"
                    )
            return

        sample_one_cols = parse_columns(args.sample_one_cols)
        sample_two_cols = parse_columns(args.sample_two_cols)
        multi_samples = [
            parse_columns(group)
            for group in (args.sample_cols or "").split(";")
            if group
        ]
        n_input_cols = len(data_lines[0].split("\t"))
        if args.output_mode == "results_only":
            fields = list(labels)
        else:
            fields = (
                list(header_fields)
                if header_fields is not None
                else ["column_%d" % i for i in range(1, n_input_cols + 1)]
            ) + labels
        if pvalue_index is not None:
            fields.append("pvalue_adjusted")
        header_row(fields)

        buffered = []
        for line_no, line in enumerate(data_lines, start=1):
            cols = line.split("\t")
            try:
                sample_one = apply_nan_policy(
                    [parse_float(cols[c - 1]) for c in sample_one_cols],
                    args.nan_policy,
                )
                sample_two = (
                    apply_nan_policy(
                        [parse_float(cols[c - 1]) for c in sample_two_cols],
                        args.nan_policy,
                    )
                    if sample_two_cols
                    else []
                )
                b_samples = [
                    apply_nan_policy(
                        [parse_float(cols[c - 1]) for c in group], args.nan_policy
                    )
                    for group in multi_samples
                ]
                results = compute(test_id, args, sample_one, sample_two, b_samples)
            except Exception:
                raise ValueError(
                    "failed to compute %s for row %d" % (test_id, line_no)
                )
            if args.output_mode == "results_only":
                out_fields = result_fields(results)
            else:
                out_fields = cols + result_fields(results)
            if pvalue_index is None:
                fout.write("\t".join(out_fields) + "\n")
            else:
                buffered.append((out_fields, results[pvalue_index]))
        if pvalue_index is not None:
            if correction == "bonferroni":
                adjusted = bonferroni([p for _, p in buffered])
            else:
                adjusted = benjamini_hochberg([p for _, p in buffered])
            for (out_fields, _), adj in zip(buffered, adjusted):
                fout.write("\t".join(out_fields + [stringify(adj)]) + "\n")


if __name__ == "__main__":
    main()
