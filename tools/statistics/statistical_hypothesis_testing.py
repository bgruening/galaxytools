#!/usr/bin/env python

"""

"""
import sys
import argparse
from scipy import stats


def columns_to_values(args, line):
    # here you go over every list
    samples = []
    for list in args:
        cols = line.split("\t")
        sample_list = []
        for row in list:
            sample_list.append(cols[row - 1])
        samples.append(map(int, sample_list))
    return samples


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--infile", required=True, help="Tabular file.")
    parser.add_argument(
        "-o", "--outfile", required=True, help="Path to the output file."
    )
    parser.add_argument("--sample_one_cols", help="Input format, like smi, sdf, inchi")
    parser.add_argument("--sample_two_cols", help="Input format, like smi, sdf, inchi")
    parser.add_argument(
        "--sample_cols",
        help="Input format, like smi, sdf, inchi,separate arrays using ;",
    )
    parser.add_argument("--test_id", help="statistical test method")
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
        help="If set perform a standard independent 2 sample test that assumes equal population variances. If not set, perform Welch's t-test, which does not assume equal population variance.",
    )
    parser.add_argument(
        "--reta",
        action="store_true",
        default=False,
        help="Whether or not to return the internally computed a values.",
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
        help="If True, if there are extra points a warning is raised saying how many of those points there are",
    )
    parser.add_argument(
        "--initial_lexsort",
        action="store_true",
        default="False",
        help="Whether to use lexsort or quicksort as the sorting method for the initial sort of the inputs.",
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
        help="Axis can equal None (ravel array first), or an integer (the axis over which to operate on a and b)",
    )
    parser.add_argument(
        "--n",
        type=int,
        default=0,
        help="the number of trials. This is ignored if x gives both the number of successes and failures",
    )
    parser.add_argument(
        "--b", type=int, default=0, help="The number of bins to use for the histogram"
    )
    parser.add_argument(
        "--N", type=int, default=0, help="Score that is compared to the elements in a."
    )
    parser.add_argument(
        "--ddof", type=int, default=0, help="Degrees of freedom correction"
    )
    parser.add_argument(
        "--score",
        type=int,
        default=0,
        help="Score that is compared to the elements in a.",
    )
    parser.add_argument("--m", type=float, default=0.0, help="limits")
    parser.add_argument("--mf", type=float, default=2.0, help="lower limit")
    parser.add_argument("--nf", type=float, default=99.9, help="higher_limit")
    parser.add_argument(
        "--p",
        type=float,
        default=0.5,
        help="The hypothesized probability of success. 0 <= p <= 1. The default value is p = 0.5",
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
        default=0,
        help="If lmbda is not None, do the transformation for that value.If lmbda is None, find the lambda that maximizes the log-likelihood function and return it as the second output argument.",
    )
    parser.add_argument(
        "--base",
        type=float,
        default=1.6,
        help="The logarithmic base to use, defaults to e",
    )
    parser.add_argument("--dtype", help="dtype")
    parser.add_argument("--med", help="med")
    parser.add_argument("--cdf", help="cdf")
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

    args = parser.parse_args()
    infile = args.infile
    outfile = open(args.outfile, "w+")
    test_id = args.test_id
    nf = args.nf
    mf = args.mf
    imbda = args.imbda
    inclusive1 = args.inclusive1
    inclusive2 = args.inclusive2
    sample0 = 0
    sample1 = 0
    sample2 = 0
    if args.sample_cols != None:
        sample0 = 1
        barlett_samples = []
        for sample in args.sample_cols.split(";"):
            barlett_samples.append(map(int, sample.split(",")))
    if args.sample_one_cols != None:
        sample1 = 1
        sample_one_cols = args.sample_one_cols.split(",")
    if args.sample_two_cols != None:
        sample_two_cols = args.sample_two_cols.split(",")
        sample2 = 1
    for line in open(infile):
        sample_one = []
        sample_two = []
        cols = line.strip().split("\t")
        if sample0 == 1:
            b_samples = columns_to_values(barlett_samples, line)
        if sample1 == 1:
            for index in sample_one_cols:
                sample_one.append(cols[int(index) - 1])
        if sample2 == 1:
            for index in sample_two_cols:
                sample_two.append(cols[int(index) - 1])
        if test_id.strip() == "describe":
            size, min_max, mean, uv, bs, bk = stats.describe(map(float, sample_one))
            cols.append(size)
            cols.append(min_max)
            cols.append(mean)
            cols.append(uv)
            cols.append(bs)
            cols.append(bk)
        elif test_id.strip() == "mode":
            vals, counts = stats.mode(map(float, sample_one))
            cols.append(vals)
            cols.append(counts)
        elif test_id.strip() == "nanmean":
            m = stats.nanmean(map(float, sample_one))
            cols.append(m)
        elif test_id.strip() == "nanmedian":
            m = stats.nanmedian(map(float, sample_one))
            cols.append(m)
        elif test_id.strip() == "kurtosistest":
            z_value, p_value = stats.kurtosistest(map(float, sample_one))
            cols.append(z_value)
            cols.append(p_value)
        elif test_id.strip() == "variation":
            ra = stats.variation(map(float, sample_one))
            cols.append(ra)
        elif test_id.strip() == "itemfreq":
            freq = stats.itemfreq(map(float, sample_one))
            for list in freq:
                elements = ",".join(map(str, list))
                cols.append(elements)
        elif test_id.strip() == "nanmedian":
            m = stats.nanmedian(map(float, sample_one))
            cols.append(m)
        elif test_id.strip() == "variation":
            ra = stats.variation(map(float, sample_one))
            cols.append(ra)
        elif test_id.strip() == "boxcox_llf":
            IIf = stats.boxcox_llf(imbda, map(float, sample_one))
            cols.append(IIf)
        elif test_id.strip() == "tiecorrect":
            fa = stats.tiecorrect(map(float, sample_one))
            cols.append(fa)
        elif test_id.strip() == "rankdata":
            r = stats.rankdata(map(float, sample_one), method=args.md)
            cols.append(r)
        elif test_id.strip() == "nanstd":
            s = stats.nanstd(map(float, sample_one), bias=args.bias)
            cols.append(s)
        elif test_id.strip() == "anderson":
            A2, critical, sig = stats.anderson(map(float, sample_one), dist=args.dist)
            cols.append(A2)
            for list in critical:
                cols.append(list)
            cols.append(",")
            for list in sig:
                cols.append(list)
        elif test_id.strip() == "binom_test":
            p_value = stats.binom_test(map(float, sample_one), n=args.n, p=args.p)
            cols.append(p_value)
        elif test_id.strip() == "gmean":
            gm = stats.gmean(map(float, sample_one), dtype=args.dtype)
            cols.append(gm)
        elif test_id.strip() == "hmean":
            hm = stats.hmean(map(float, sample_one), dtype=args.dtype)
            cols.append(hm)
        elif test_id.strip() == "kurtosis":
            k = stats.kurtosis(
                map(float, sample_one),
                axis=args.axis,
                fisher=args.fisher,
                bias=args.bias,
            )
            cols.append(k)
        elif test_id.strip() == "moment":
            n_moment = stats.moment(map(float, sample_one), n=args.n)
            cols.append(n_moment)
        elif test_id.strip() == "normaltest":
            k2, p_value = stats.normaltest(map(float, sample_one))
            cols.append(k2)
            cols.append(p_value)
        elif test_id.strip() == "skew":
            skewness = stats.skew(map(float, sample_one), bias=args.bias)
            cols.append(skewness)
        elif test_id.strip() == "skewtest":
            z_value, p_value = stats.skewtest(map(float, sample_one))
            cols.append(z_value)
            cols.append(p_value)
        elif test_id.strip() == "sem":
            s = stats.sem(map(float, sample_one), ddof=args.ddof)
            cols.append(s)
        elif test_id.strip() == "zscore":
            z = stats.zscore(map(float, sample_one), ddof=args.ddof)
            for list in z:
                cols.append(list)
        elif test_id.strip() == "signaltonoise":
            s2n = stats.signaltonoise(map(float, sample_one), ddof=args.ddof)
            cols.append(s2n)
        elif test_id.strip() == "percentileofscore":
            p = stats.percentileofscore(
                map(float, sample_one), score=args.score, kind=args.kind
            )
            cols.append(p)
        elif test_id.strip() == "bayes_mvs":
            c_mean, c_var, c_std = stats.bayes_mvs(
                map(float, sample_one), alpha=args.alpha
            )
            cols.append(c_mean)
            cols.append(c_var)
            cols.append(c_std)
        elif test_id.strip() == "sigmaclip":
            c, c_low, c_up = stats.sigmaclip(
                map(float, sample_one), low=args.m, high=args.n
            )
            cols.append(c)
            cols.append(c_low)
            cols.append(c_up)
        elif test_id.strip() == "kstest":
            d, p_value = stats.kstest(
                map(float, sample_one),
                cdf=args.cdf,
                N=args.N,
                alternative=args.alternative,
                mode=args.mode,
            )
            cols.append(d)
            cols.append(p_value)
        elif test_id.strip() == "chi2_contingency":
            chi2, p, dof, ex = stats.chi2_contingency(
                map(float, sample_one), correction=args.correction, lambda_=args.lambda_
            )
            cols.append(chi2)
            cols.append(p)
            cols.append(dof)
            cols.append(ex)
        elif test_id.strip() == "tmean":
            if nf is 0 and mf is 0:
                mean = stats.tmean(map(float, sample_one))
            else:
                mean = stats.tmean(
                    map(float, sample_one), (mf, nf), (inclusive1, inclusive2)
                )
            cols.append(mean)
        elif test_id.strip() == "tmin":
            if mf is 0:
                min = stats.tmin(map(float, sample_one))
            else:
                min = stats.tmin(
                    map(float, sample_one), lowerlimit=mf, inclusive=args.inclusive
                )
            cols.append(min)
        elif test_id.strip() == "tmax":
            if nf is 0:
                max = stats.tmax(map(float, sample_one))
            else:
                max = stats.tmax(
                    map(float, sample_one), upperlimit=nf, inclusive=args.inclusive
                )
            cols.append(max)
        elif test_id.strip() == "tvar":
            if nf is 0 and mf is 0:
                var = stats.tvar(map(float, sample_one))
            else:
                var = stats.tvar(
                    map(float, sample_one), (mf, nf), (inclusive1, inclusive2)
                )
            cols.append(var)
        elif test_id.strip() == "tstd":
            if nf is 0 and mf is 0:
                std = stats.tstd(map(float, sample_one))
            else:
                std = stats.tstd(
                    map(float, sample_one), (mf, nf), (inclusive1, inclusive2)
                )
            cols.append(std)
        elif test_id.strip() == "tsem":
            if nf is 0 and mf is 0:
                s = stats.tsem(map(float, sample_one))
            else:
                s = stats.tsem(
                    map(float, sample_one), (mf, nf), (inclusive1, inclusive2)
                )
            cols.append(s)
        elif test_id.strip() == "scoreatpercentile":
            if nf is 0 and mf is 0:
                s = stats.scoreatpercentile(
                    map(float, sample_one),
                    map(float, sample_two),
                    interpolation_method=args.interpolation,
                )
            else:
                s = stats.scoreatpercentile(
                    map(float, sample_one),
                    map(float, sample_two),
                    (mf, nf),
                    interpolation_method=args.interpolation,
                )
            for list in s:
                cols.append(list)
        elif test_id.strip() == "relfreq":
            if nf is 0 and mf is 0:
                rel, low_range, binsize, ex = stats.relfreq(
                    map(float, sample_one), args.b
                )
            else:
                rel, low_range, binsize, ex = stats.relfreq(
                    map(float, sample_one), args.b, (mf, nf)
                )
            for list in rel:
                cols.append(list)
            cols.append(low_range)
            cols.append(binsize)
            cols.append(ex)
        elif test_id.strip() == "binned_statistic":
            if nf is 0 and mf is 0:
                st, b_edge, b_n = stats.binned_statistic(
                    map(float, sample_one),
                    map(float, sample_two),
                    statistic=args.statistic,
                    bins=args.b,
                )
            else:
                st, b_edge, b_n = stats.binned_statistic(
                    map(float, sample_one),
                    map(float, sample_two),
                    statistic=args.statistic,
                    bins=args.b,
                    range=(mf, nf),
                )
            cols.append(st)
            cols.append(b_edge)
            cols.append(b_n)
        elif test_id.strip() == "threshold":
            if nf is 0 and mf is 0:
                o = stats.threshold(map(float, sample_one), newval=args.new)
            else:
                o = stats.threshold(map(float, sample_one), mf, nf, newval=args.new)
            for list in o:
                cols.append(list)
        elif test_id.strip() == "trimboth":
            o = stats.trimboth(
                map(float, sample_one), proportiontocut=args.proportiontocut
            )
            for list in o:
                cols.append(list)
        elif test_id.strip() == "trim1":
            t1 = stats.trim1(
                map(float, sample_one),
                proportiontocut=args.proportiontocut,
                tail=args.tail,
            )
            for list in t1:
                cols.append(list)
        elif test_id.strip() == "histogram":
            if nf is 0 and mf is 0:
                hi, low_range, binsize, ex = stats.histogram(
                    map(float, sample_one), args.b
                )
            else:
                hi, low_range, binsize, ex = stats.histogram(
                    map(float, sample_one), args.b, (mf, nf)
                )
            cols.append(hi)
            cols.append(low_range)
            cols.append(binsize)
            cols.append(ex)
        elif test_id.strip() == "cumfreq":
            if nf is 0 and mf is 0:
                cum, low_range, binsize, ex = stats.cumfreq(
                    map(float, sample_one), args.b
                )
            else:
                cum, low_range, binsize, ex = stats.cumfreq(
                    map(float, sample_one), args.b, (mf, nf)
                )
            cols.append(cum)
            cols.append(low_range)
            cols.append(binsize)
            cols.append(ex)
        elif test_id.strip() == "boxcox_normmax":
            if nf is 0 and mf is 0:
                ma = stats.boxcox_normmax(map(float, sample_one))
            else:
                ma = stats.boxcox_normmax(
                    map(float, sample_one), (mf, nf), method=args.method
                )
            cols.append(ma)
        elif test_id.strip() == "boxcox":
            if imbda is 0:
                box, ma, ci = stats.boxcox(map(float, sample_one), alpha=args.alpha)
                cols.append(box)
                cols.append(ma)
                cols.append(ci)
            else:
                box = stats.boxcox(map(float, sample_one), imbda, alpha=args.alpha)
                cols.append(box)
        elif test_id.strip() == "histogram2":
            h2 = stats.histogram2(map(float, sample_one), map(float, sample_two))
            for list in h2:
                cols.append(list)
        elif test_id.strip() == "ranksums":
            z_statistic, p_value = stats.ranksums(
                map(float, sample_one), map(float, sample_two)
            )
            cols.append(z_statistic)
            cols.append(p_value)
        elif test_id.strip() == "ttest_1samp":
            t, prob = stats.ttest_1samp(map(float, sample_one), map(float, sample_two))
            for list in t:
                cols.append(list)
            for list in prob:
                cols.append(list)
        elif test_id.strip() == "ansari":
            AB, p_value = stats.ansari(map(float, sample_one), map(float, sample_two))
            cols.append(AB)
            cols.append(p_value)
        elif test_id.strip() == "linregress":
            slope, intercept, r_value, p_value, stderr = stats.linregress(
                map(float, sample_one), map(float, sample_two)
            )
            cols.append(slope)
            cols.append(intercept)
            cols.append(r_value)
            cols.append(p_value)
            cols.append(stderr)
        elif test_id.strip() == "pearsonr":
            cor, p_value = stats.pearsonr(
                map(float, sample_one), map(float, sample_two)
            )
            cols.append(cor)
            cols.append(p_value)
        elif test_id.strip() == "pointbiserialr":
            r, p_value = stats.pointbiserialr(
                map(float, sample_one), map(float, sample_two)
            )
            cols.append(r)
            cols.append(p_value)
        elif test_id.strip() == "ks_2samp":
            d, p_value = stats.ks_2samp(map(float, sample_one), map(float, sample_two))
            cols.append(d)
            cols.append(p_value)
        elif test_id.strip() == "mannwhitneyu":
            mw_stats_u, p_value = stats.mannwhitneyu(
                map(float, sample_one),
                map(float, sample_two),
                use_continuity=args.mwu_use_continuity,
            )
            cols.append(mw_stats_u)
            cols.append(p_value)
        elif test_id.strip() == "zmap":
            z = stats.zmap(
                map(float, sample_one), map(float, sample_two), ddof=args.ddof
            )
            for list in z:
                cols.append(list)
        elif test_id.strip() == "ttest_ind":
            mw_stats_u, p_value = stats.ttest_ind(
                map(float, sample_one), map(float, sample_two), equal_var=args.equal_var
            )
            cols.append(mw_stats_u)
            cols.append(p_value)
        elif test_id.strip() == "ttest_rel":
            t, prob = stats.ttest_rel(
                map(float, sample_one), map(float, sample_two), axis=args.axis
            )
            cols.append(t)
            cols.append(prob)
        elif test_id.strip() == "mood":
            z, p_value = stats.mood(
                map(float, sample_one), map(float, sample_two), axis=args.axis
            )
            cols.append(z)
            cols.append(p_value)
        elif test_id.strip() == "shapiro":
            W, p_value, a = stats.shapiro(
                map(float, sample_one), map(float, sample_two), args.reta
            )
            cols.append(W)
            cols.append(p_value)
            for list in a:
                cols.append(list)
        elif test_id.strip() == "kendalltau":
            k, p_value = stats.kendalltau(
                map(float, sample_one),
                map(float, sample_two),
                initial_lexsort=args.initial_lexsort,
            )
            cols.append(k)
            cols.append(p_value)
        elif test_id.strip() == "entropy":
            s = stats.entropy(
                map(float, sample_one), map(float, sample_two), base=args.base
            )
            cols.append(s)
        elif test_id.strip() == "spearmanr":
            if sample2 == 1:
                rho, p_value = stats.spearmanr(
                    map(float, sample_one), map(float, sample_two)
                )
            else:
                rho, p_value = stats.spearmanr(map(float, sample_one))
            cols.append(rho)
            cols.append(p_value)
        elif test_id.strip() == "wilcoxon":
            if sample2 == 1:
                T, p_value = stats.wilcoxon(
                    map(float, sample_one),
                    map(float, sample_two),
                    zero_method=args.zero_method,
                    correction=args.correction,
                )
            else:
                T, p_value = stats.wilcoxon(
                    map(float, sample_one),
                    zero_method=args.zero_method,
                    correction=args.correction,
                )
            cols.append(T)
            cols.append(p_value)
        elif test_id.strip() == "chisquare":
            if sample2 == 1:
                rho, p_value = stats.chisquare(
                    map(float, sample_one), map(float, sample_two), ddof=args.ddof
                )
            else:
                rho, p_value = stats.chisquare(map(float, sample_one), ddof=args.ddof)
            cols.append(rho)
            cols.append(p_value)
        elif test_id.strip() == "power_divergence":
            if sample2 == 1:
                stat, p_value = stats.power_divergence(
                    map(float, sample_one),
                    map(float, sample_two),
                    ddof=args.ddof,
                    lambda_=args.lambda_,
                )
            else:
                stat, p_value = stats.power_divergence(
                    map(float, sample_one), ddof=args.ddof, lambda_=args.lambda_
                )
            cols.append(stat)
            cols.append(p_value)
        elif test_id.strip() == "theilslopes":
            if sample2 == 1:
                mpe, met, lo, up = stats.theilslopes(
                    map(float, sample_one), map(float, sample_two), alpha=args.alpha
                )
            else:
                mpe, met, lo, up = stats.theilslopes(
                    map(float, sample_one), alpha=args.alpha
                )
            cols.append(mpe)
            cols.append(met)
            cols.append(lo)
            cols.append(up)
        elif test_id.strip() == "combine_pvalues":
            if sample2 == 1:
                stat, p_value = stats.combine_pvalues(
                    map(float, sample_one),
                    method=args.med,
                    weights=map(float, sample_two),
                )
            else:
                stat, p_value = stats.combine_pvalues(
                    map(float, sample_one), method=args.med
                )
            cols.append(stat)
            cols.append(p_value)
        elif test_id.strip() == "obrientransform":
            ob = stats.obrientransform(*b_samples)
            for list in ob:
                elements = ",".join(map(str, list))
                cols.append(elements)
        elif test_id.strip() == "f_oneway":
            f_value, p_value = stats.f_oneway(*b_samples)
            cols.append(f_value)
            cols.append(p_value)
        elif test_id.strip() == "kruskal":
            h, p_value = stats.kruskal(*b_samples)
            cols.append(h)
            cols.append(p_value)
        elif test_id.strip() == "friedmanchisquare":
            fr, p_value = stats.friedmanchisquare(*b_samples)
            cols.append(fr)
            cols.append(p_value)
        elif test_id.strip() == "fligner":
            xsq, p_value = stats.fligner(
                center=args.center, proportiontocut=args.proportiontocut, *b_samples
            )
            cols.append(xsq)
            cols.append(p_value)
        elif test_id.strip() == "bartlett":
            T, p_value = stats.bartlett(*b_samples)
            cols.append(T)
            cols.append(p_value)
        elif test_id.strip() == "levene":
            w, p_value = stats.levene(
                center=args.center, proportiontocut=args.proportiontocut, *b_samples
            )
            cols.append(w)
            cols.append(p_value)
        elif test_id.strip() == "median_test":
            stat, p_value, m, table = stats.median_test(
                ties=args.ties,
                correction=args.correction,
                lambda_=args.lambda_,
                *b_samples
            )
            cols.append(stat)
            cols.append(p_value)
            cols.append(m)
            cols.append(table)
            for list in table:
                elements = ",".join(map(str, list))
                cols.append(elements)
        outfile.write("%s\n" % "\t".join(map(str, cols)))
    outfile.close()


if __name__ == "__main__":
    main()
