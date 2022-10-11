#!/usr/bin/env python3

import argparse
import os
import subprocess


"""
DEWSeq wrapper script dependencies:
conda install -c bioconda bioconductor-dewseq
conda install -c conda-forge r-rmarkdown
conda install -c bioconda bioconductor-biocstyle
conda install -c conda-forge r-tidyverse
conda install -c bioconda r-ggrepel
conda install -c bioconda bioconductor-ihw

Wrapper for DEWSeq R markdown file:
https://github.com/EMBL-Hentze-group/DEWSeq_analysis_helpers/tree/master/Parametrized_Rmd


Test runs
=========

# This reports 150 significant regions.
python dewseq_wrapper.py --annot test-data/windows.exp.txt --matrix test-data/Rbp_count_matrix.exp.txt --info test-data/sample_info.exp.txt --out test_out --ds-pvc 0.5

# Wherease, with LRT, DEWSeq reports just one siginificant region.
python dewseq_wrapper.py --annot test-data/windows.exp.txt --matrix test-data/Rbp_count_matrix.exp.txt --info test-data/sample_info.exp.txt --out test2_out --ds-pvc 0.5 --ds-use-lrt

"""

################################################################################


def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Wrapping DEWSeq R markdown script, to call peak regions on the CLIP-seq data,
    preprocessed by htseq-clip.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="dewseq_wrapper.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--annot",
                   dest="in_annot",
                   type=str,
                   metavar='str',
                   required=True,
                   help="DEWseq annotation file, i.e. windows mapped to IDs table file (output of htseq-clip mapToId)")
    p.add_argument("--matrix",
                   dest="in_matrix",
                   type=str,
                   metavar='str',
                   required=True,
                   help="DEWseq count matrix file (output of htseq-clip createMatrix)")
    p.add_argument("--info",
                   dest="in_info",
                   type=str,
                   metavar='str',
                   required=True,
                   help="DEWseq sample information file (output of htsc_create_count_table.py / htseq-clip Create count table Galaxy wrapper)")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Results output folder")
    p.add_argument("--ds-ms",
                   dest="ds_ms",
                   type=int,
                   metavar='int',
                   default=2,
                   help="DEWSeq min_sample parameter. Keep only the windows with at least min_sample number of samples with crosslink site count > min_count (default: 2)")
    p.add_argument("--ds-mc",
                   dest="ds_mc",
                   type=int,
                   metavar='int',
                   default=2,
                   help="DEWSeq min_count parameter. Minimum crosslink site per window per sample (default: 2)")
    p.add_argument("--ds-pvc",
                   dest="ds_pvc",
                   type=float,
                   metavar='float',
                   default=0.1,
                   help="DEWSeq p_value_cutoff parameter. p adjusted value threshold for significant windows (default: 0.1)")
    p.add_argument("--ds-lfcc",
                   dest="ds_lfcc",
                   type=float,
                   metavar='float',
                   default=1,
                   help="DEWSeq lfc_cutoff parameter. Log2 fold change threshold for significant windows (default: 1)")
    p.add_argument("--ds-oc",
                   dest="ds_use_oc",
                   default=False,
                   action="store_true",
                   help="DEWSeq overlap_correction parameter. By default FALSE, i.e., do not adjust p-value for overlapping windows. If TRUE use Bonferroni family wise error rate correction on overlapping sliding windows (default: false)")
    p.add_argument("--ds--disable-ihw",
                   dest="ds_disable_ihw",
                   default=False,
                   action="store_true",
                   help="Disable DEWSeq IHW parameter. By default, use IHW for multiple tesing correction instead of default BH (Benjamini Hochberg). NOTE: We recommend using IHW instead of default BH for FDR correction")
    p.add_argument("--ds--disable-df",
                   dest="ds_disable_df",
                   default=False,
                   action="store_true",
                   help="Disable DEWSeq decide_fit parameter. By default, decide on dispersion estimation fit type local or parametric. If disabled, Use parametric fit. NOTE: decide_fit=TRUE will fit data using both parametric and local fit types and will choose the best fit of the two (see DESeq2 vignette for details). Typically, this should give better results compared to using the default fit type parametric. But, keep in mind that this will also increase the total run time")
    p.add_argument("--ds-use-lrt",
                   dest="ds_use_lrt",
                   default=False,
                   action="store_true",
                   help="DEWSeq LRT parameter. Use LRT if the given value is TRUE (see DESeq2 vignette for details). By default, DEWSeq uses Wald test.  NOTE: In our experience, LRT is more accurate than Wald test. But, keep in mind that LRT is a stringent test in comparison to Wald. So if your protein of interest is a very active binder, run the analysis with LRT=TRUE, otherwise use it with caution as you may end up with no significant windows or regions in your final output")
    p.add_argument("--ds-id",
                   dest="ds_id",
                   type=str,
                   default="RBP",
                   metavar='str',
                   help="DEWSeq dataset ID for output report (default: RBP)")
    p.add_argument("--ds-markdown",
                   dest="ds_markdown",
                   type=str,
                   default="analyseStudy.Rmd",
                   metavar='str',
                   help="Provide path to DEWSeq markdown file. By default assumed to be in working directory")
    return p


################################################################################

def count_file_rows(in_file,
                    nr_cols=False):
    """
    Count number of file rows. If nr_cols set, demand certain (nr_cols) number
    of columns (separated by tab), in order for row to be counted.

    """
    c = 0
    with open(in_file) as f:
        for line in f:
            cols = line.strip().split("\t")
            if nr_cols:
                if len(cols) == nr_cols:
                    c += 1
            else:
                c += 1
    f.closed
    return c


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_annot), "--annot file \"%s\" not found" % (args.in_annot)
    assert os.path.exists(args.in_matrix), "--matrix file \"%s\" not found" % (args.in_matrix)
    assert os.path.exists(args.in_info), "--info file \"%s\" not found" % (args.in_info)
    assert os.path.exists(args.ds_markdown), "--ds-markdown file \"%s\" not found" % (args.ds_markdown)

    # Input files.
    annot_in = os.path.abspath(args.in_annot)
    matrix_in = os.path.abspath(args.in_matrix)
    info_in = os.path.abspath(args.in_info)
    md_in = os.path.abspath(args.ds_markdown)

    # Output folder.
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    # Output files.
    abs_path_out = os.path.abspath(args.out_folder)

    html_out = abs_path_out + "/report.html"
    win_csv_out = abs_path_out + "/windows.csv"
    sig_reg_csv_out = abs_path_out + "/significant_regions.csv"
    sig_win_reg_bed_out = abs_path_out + "/significant_windows_and_regions.bed"
    sig_reg_bed_out = abs_path_out + "/significant_regions.bed"

    # Delete existing files (as if no peaks found old files would be reported).
    if os.path.exists(html_out):
        os.remove(html_out)
    if os.path.exists(win_csv_out):
        os.remove(win_csv_out)
    if os.path.exists(sig_reg_csv_out):
        os.remove(sig_reg_csv_out)
    if os.path.exists(sig_win_reg_bed_out):
        os.remove(sig_win_reg_bed_out)
    if os.path.exists(sig_reg_bed_out):
        os.remove(sig_reg_bed_out)

    """
    Run DEWSeq R markdown file.

    """

    md_ihw = "TRUE"
    md_df = "TRUE"
    md_lrt = "FALSE"
    md_olc = "FALSE"
    if args.ds_disable_ihw:
        md_ihw = "FALSE"
    if args.ds_disable_df:
        md_df = "FALSE"
    if args.ds_use_lrt:
        md_lrt = "TRUE"
    if args.ds_use_oc:
        md_olc = "TRUE"

    md_cmd = "Rscript -e 'rmarkdown::render("
    md_cmd += 'input = "%s", ' % (md_in)
    md_cmd += 'output_file = "%s", ' % (html_out)
    md_cmd += 'params = list(protein = "%s", ' % (args.ds_id)
    md_cmd += 'sampleinfo_file = "%s", ' % (info_in)
    md_cmd += 'countmatrix_file = "%s", ' % (matrix_in)
    md_cmd += 'annotation_file = "%s", ' % (annot_in)
    md_cmd += 'output_windows_file = "%s", ' % (win_csv_out)
    md_cmd += 'output_regions_file = "%s", ' % (sig_reg_csv_out)
    md_cmd += 'output_bed_file = "%s", ' % (sig_win_reg_bed_out)
    md_cmd += 'min_count = %i, ' % (args.ds_mc)
    md_cmd += 'min_sample = %i, ' % (args.ds_ms)
    md_cmd += 'p_value_cutoff = %s, ' % (str(args.ds_pvc))
    md_cmd += 'lfc_cutoff = %s, ' % (str(args.ds_lfcc))
    md_cmd += 'overlap_correction = %s, ' % (md_olc)
    md_cmd += 'IHW = %s, ' % (md_ihw)
    md_cmd += 'decide_fit = %s, ' % (md_df)
    md_cmd += 'LRT = %s))' % (md_lrt)
    md_cmd += "'"

    print("Running DEWSeq R markdown file ... ")
    print(md_cmd)
    output = subprocess.getoutput(md_cmd)
    print(output)

    assert os.path.exists(html_out), "output file \"%s\" not found" % (html_out)
    assert os.path.exists(win_csv_out), "output file \"%s\" not found" % (win_csv_out)

    print("")
    if not os.path.exists(sig_reg_csv_out):
        print("WARNING: no significant regions found! (missing \"%s\" file)" % (sig_reg_csv_out))
    else:
        assert os.path.exists(sig_win_reg_bed_out), "output file \"%s\" not found" % (sig_win_reg_bed_out)
        c_sig_reg = count_file_rows(sig_reg_csv_out) - 1
        print("# significant regions:  %i" % (c_sig_reg))
        # Save contiguous BED regions in separate file.
        OUTBED = open(sig_reg_bed_out, "w")
        with open(sig_reg_csv_out) as f:
            for line in f:
                row = line.strip()
                cols = line.strip().split("\t")
                if cols[1] == "region_begin":
                    continue
                chr_id = cols[0]
                reg_s = cols[1]
                reg_e = cols[2]
                pol = cols[3]
                win_in_reg = int(cols[4])
                padj_mean = cols[7]
                logfc_mean = cols[10]
                reg_id = cols[12]
                new_reg_id = reg_id
                # For regions consisting of > 1 window.
                if win_in_reg > 1:
                    new_reg_id = reg_id + "@region"
                # Print out BED with additional columns (padj, logfc).
                OUTBED.write("%s\t%s\t%s\t%s\t0\t%s\t%s\t%s\n" % (chr_id, reg_s, reg_e, new_reg_id, pol, padj_mean, logfc_mean))

        f.closed
        OUTBED.close()

    """
    Report.

    """

    print("")
    print("OUTPUT FILES")
    print("============")
    print("HTML report:\n%s" % (html_out))
    print("Windows CSV:\n%s" % (win_csv_out))
    if os.path.exists(sig_reg_csv_out):
        print("Significant regions CSV:\n%s" % (sig_reg_csv_out))
    if os.path.exists(sig_win_reg_bed_out):
        print("Significant windows + regions BED:\n%s" % (sig_win_reg_bed_out))
    if os.path.exists(sig_reg_bed_out):
        print("Significant regions BED:\n%s" % (sig_reg_bed_out))
    print("")
