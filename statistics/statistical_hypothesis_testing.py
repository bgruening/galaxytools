#!/usr/bin/env python

"""
    
"""
import sys
import argparse
from scipy.stats import ranksums
from scipy.stats import mannwhitneyu
from scipy.stats import ttest_ind



def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Tabular file.')
    parser.add_argument('-o', '--outfile', required=True, help='Path to the output file.')
    parser.add_argument("--sample_one_cols", help="Input format, like smi, sdf, inchi")
    parser.add_argument("--sample_two_cols", help="Input format, like smi, sdf, inchi")
    parser.add_argument("--test_id", help="statistical test method")

    parser.add_argument("--mwu_use_continuity", action="store_true", default="False",
                    help="Whether a continuity correction (1/2.) should be taken into account.")
    parser.add_argument("--equal_var", action="store_true", default="False",
                    help="If set perform a standard independent 2 sample test that assumes equal population variances. If not set, perform Welch's t-test, which does not assume equal population variance.")


    args = parser.parse_args()

    infile = args.infile
    outfile = open(args.outfile, 'w+')
    test_id = args.test_id
    sample_one_cols = args.sample_one_cols.split(',')
    sample_two_cols = args.sample_two_cols.split(',')

    for line in open( infile ):
        cols = line.strip().split('\t')
        sample_one = []
        sample_two = []
        for index in sample_one_cols:
            sample_one.append( cols[ int(index) -1 ] )
        for index in sample_two_cols:
            sample_two.append( cols[ int(index) -1 ] )

        if test_id.strip() == 'ranksums':
            z_statistic, p_value = ranksums( map(float, sample_one), map(float, sample_two) )
            cols.append(z_statistic)
            cols.append(p_value)
        elif test_id.strip() == 'mannwhitneyu':
            mw_stats_u, p_value = mannwhitneyu( map(float, sample_one), map(float, sample_two), use_continuity=args.mwu_use_continuity )
            cols.append( mw_stats_u )
            cols.append( p_value )
        elif test_id.strip() == 'ttest_ind':
            mw_stats_u, p_value = ttest_ind( map(float, sample_one), map(float, sample_two), equal_var=args.equal_var )
            cols.append( mw_stats_u )
            cols.append( p_value )

        outfile.write( '%s\n' % '\t'.join( map(str, cols) ) )
    outfile.close()


if __name__ == '__main__':
    main()
