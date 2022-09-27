#!/usr/bin/env python3

import argparse
import os
import subprocess

import pysam

"""

Install deseq-clip
==================

conda install -c bioconda pysam
conda install -c bioconda htseq
pip install htseq-clip

Or directly by:
conda install -c bioconda htseq-clip


Test call
=========

python htsc_create_count_table.py --data-id Rbp --win-bed test-data/windows.exp.bed --out test_create_count_table_out --exp-bams test-data/Rbp_exp_rep1.Synechocystis_pSYSM.bam test-data/Rbp_exp_rep2.Synechocystis_pSYSM.bam --ctr-bams test-data/Rbp_ctrl_rep1.Synechocystis_pSYSM.bam --no-zipper

Compare:
diff test-data/Rbp_count_matrix.exp.txt test_create_count_table_out/count_matrix.txt
diff test-data/sample_info.exp.txt test_create_count_table_out/sample_info.txt

This corresponds to:
htseq-clip extract -i test-data/Rbp_exp_rep1.Synechocystis_pSYSM.bam -o test_create_count_table_out/Rbp_exp_rep1.Synechocystis_pSYSM.bed -e 1 -s m -g 0 -q 10 -c 4 -m 0 -x 500 -l 10000
htseq-clip extract -i test-data/Rbp_exp_rep2.Synechocystis_pSYSM.bam -o test_create_count_table_out/Rbp_exp_rep2.Synechocystis_pSYSM.bed -e 1 -s m -g 0 -q 10 -c 4 -m 0 -x 500 -l 10000
htseq-clip extract -i test-data/Rbp_ctrl_rep1.Synechocystis_pSYSM.bam -o test_create_count_table_out/Rbp_ctrl_rep1.Synechocystis_pSYSM.bed -e 1 -s m -g 0 -q 10 -c 4 -m 0 -x 500 -l 10000
htseq-clip count -i test_create_count_table_out/Rbp_exp_rep1.Synechocystis_pSYSM.bed -a test-data/windows.exp.bed -o test_create_count_table_out/counts_Rbp/Rbp_exp_rep1.Synechocystis_pSYSM.csv
htseq-clip count -i test_create_count_table_out/Rbp_exp_rep2.Synechocystis_pSYSM.bed -a test-data/windows.exp.bed -o test_create_count_table_out/counts_Rbp/Rbp_exp_rep2.Synechocystis_pSYSM.csv
htseq-clip count -i test_create_count_table_out/Rbp_ctrl_rep1.Synechocystis_pSYSM.bed -a test-data/windows.exp.bed -o test_create_count_table_out/counts_Rbp/Rbp_ctrl_rep1.Synechocystis_pSYSM.csv
htseq-clip createMatrix -i test_create_count_table_out/counts_Rbp -o test_create_count_table_out/Rbp_count_matrix.txt -b Rbp


To get BAM content for single chromosome:
samtools view -b -h bam_all/Rbp3_uv_rep1.bam Synechocystis_pSYSM > test-data/Rbp_exp_rep1.Synechocystis_pSYSM.bam

More tests:
python htsc_create_count_table.py --hce-f test-data/chr_names.txt --data-id Rbp --win-bed test-data/windows.exp.bed --out test_create_count_table_out --exp-bams test-data/Rbp_exp_rep1.Synechocystis_pSYSM.bam test-data/Rbp_exp_rep2.Synechocystis_pSYSM.bam --ctr-bams test-data/Rbp_ctrl_rep1.Synechocystis_pSYSM.bam --no-zipper

DEWSeq input files:
test_create_count_table_out/Rbp_count_matrix.txt
test_create_count_table_out/sample_info.txt

"""


################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Use htseq-clip to extract crosslink sites from BAM files, and count
    overlaps with window regions. Finally create an R count matrix, e.g. as
    input for DEWSeq.
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="htsc_create_count_table.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--win-bed",
                   dest="in_win_bed",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Sliding windows BED annotation file created with htseq-clip createSlidingWindows (also accepts .bed.gz)")
    p.add_argument("--exp-bams",
                   dest="exp_bam_list",
                   type=str,
                   metavar='str',
                   nargs='+',
                   required=True,
                   help="List of IP BAM files (--exp-bams ip_rep1.bam ip_rep2.bam .. )")
    p.add_argument("--ctr-bams",
                   dest="ctr_bam_list",
                   type=str,
                   metavar='str',
                   nargs='+',
                   required=True,
                   help="List of control (SM input) BAM files (--ctr-bams smi_rep1.bam smi_rep2.bam .. )")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Results output folder")
    # htseq-clip extract.
    p.add_argument("--hce-e",
                   dest="hce_e",
                   type=int,
                   default=1,
                   choices=[1, 2],
                   help="htseq-clip extract -e parameter. This selects read/mate to extract crosslink sites from (default: 1)")
    p.add_argument("--hce-s",
                   dest="hce_s",
                   type=str,
                   default="m",
                   help="htseq-clip extract -s parameter. Choose crosslink site (s: start, m: middle, e: end ... ) (default: m)")
    p.add_argument("--hce-g",
                   dest="hce_g",
                   type=int,
                   metavar='int',
                   default=0,
                   help="htseq-clip extract -g parameter. Number of nucleotides to offset for crosslink sites (default: 0)")
    p.add_argument("--hce-q",
                   dest="hce_q",
                   type=int,
                   metavar='int',
                   default=10,
                   help="htseq-clip extract -q parameter. Minimum alignment quality for filtering (default: 10)")
    p.add_argument("--hce-primary",
                   dest="hce_primary",
                   default=False,
                   action="store_true",
                   help="htseq-clip extract --primary parameter. Flag to use only primary positions of multimapping reads (default: False)")
    p.add_argument("--hce-c",
                   dest="hce_c",
                   type=int,
                   metavar='int',
                   default=1,
                   help="htseq-clip extract -c parameter. Number of cores (default: 1)")
    p.add_argument("--hce-m",
                   dest="hce_m",
                   type=int,
                   metavar='int',
                   default=0,
                   help="htseq-clip extract -m parameter. Minimum read length for filtering (default: 0)")
    p.add_argument("--hce-x",
                   dest="hce_x",
                   type=int,
                   metavar='int',
                   default=500,
                   help="htseq-clip extract -x parameter. Maximum read length for filtering (default: 500)")
    p.add_argument("--hce-l",
                   dest="hce_l",
                   type=int,
                   metavar='int',
                   default=10000,
                   help="htseq-clip extract -l parameter. Maximum read interval length (default: 10000)")
    p.add_argument("--hce-f",
                   dest="hce_f",
                   type=str,
                   help="htseq-clip extract -f parameter. Extract crosslink sites only from chromosomes given in this file (one chromosome ID per line)")
    # htseq-clip count.
    p.add_argument("--hcc-unstranded",
                   dest="hcc_unstranded",
                   default=False,
                   action="store_true",
                   help="htseq-clip count --unstranded parameter. Use this flag for non strand specific crosslink site counting (default: strand-specific counting)")
    # More.
    p.add_argument("--data-id",
                   dest="dataset_id",
                   type=str,
                   default="RBP",
                   metavar='str',
                   help="Dataset ID used as prefix for naming datasets (default: RBP)")
    p.add_argument("--filter-bed",
                   dest="filter_bed",
                   type=str,
                   metavar='str',
                   help="Provide BED file to filter out BAM reads overlapping with --filter-bed regions")
    p.add_argument("--filter-mode",
                   dest="filter_mode",
                   type=int,
                   default=1,
                   choices=[1, 2],
                   help="Filter mode for --filter-bed file. 1: keep BAM reads not overlapping with --filter-bed regions. 2: Keep only BAM reads overlapping with --filter-bed (default: 1)")
    p.add_argument("--no-zipper",
                   dest="no_zipper",
                   default=False,
                   action="store_true",
                   help="Do not gzip output files (default: False)")
    p.add_argument("--keep-imf",
                   dest="keep_imf",
                   default=False,
                   action="store_true",
                   help="Keep intermediate files (filtered BAM, BED, CSV) (default: False)")
    return p


################################################################################

def bam_remove_overlap_region_reads(in_bed, in_bam, out_bam,
                                    params="",
                                    sort_bed=False):

    """
    Remove BAM reads from in_bam, based on overlap with in_bed BED regions.
    I.e. only BAM reads not overlapping with in_bed regions are stored in
    out_bam.

    Using intersectBed instead of samtools view -L, which allows -s for strand
    specific filtering.

    """
    assert os.path.exists(in_bed), "in_bed does not exist"
    assert os.path.exists(in_bam), "in_bam does not exist"

    if sort_bed:
        check_cmd = "sort -k1,1 -k2,2n " + in_bed + " | " + "intersectBed -abam " + in_bam + " -b stdin " + params + " -sorted > " + out_bam
    else:
        check_cmd = "intersectBed -abam " + in_bam + " -b " + in_bed + " " + params + " > " + out_bam
    output = subprocess.getoutput(check_cmd)
    error = False
    if output:
        error = True
    assert error is False, "intersectBed has problems with your input:\n%s\n%s" % (check_cmd, output)


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_win_bed), "--win-bed file \"%s\" not found" % (args.in_win_bed)

    for bam_file in args.exp_bam_list:
        assert os.path.exists(bam_file), "--exp-bams BAM file \"%s\" not found" % (bam_file)

    for bam_file in args.ctr_bam_list:
        assert os.path.exists(bam_file), "--ctr-bams BAM file \"%s\" not found" % (bam_file)

    # assert len(args.exp_bam_list) > 1, "# --exp-bams needs to be > 1"

    if args.filter_bed:
        assert os.path.exists(args.filter_bed), "--filter-bed file \"%s\" not found" % (args.filter_bed)

    hce_s_dic = {"s": 1, "i": 1, "d": 1, "m": 1, "e": 1}
    assert args.hce_s in hce_s_dic, "invalid --hce-s given (choose between {s,i,d,m,e})"

    # Crop dataset ID.
    data_id = args.dataset_id
    if len(args.dataset_id) > 20:
        data_id = args.dataset_id[:20]
    # Remove spaces.
    data_id = data_id.replace(" ", "_")

    # Output folders.
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)
    counts_folder = args.out_folder + "/counts_%s" % (data_id)
    if not os.path.exists(counts_folder):
        os.makedirs(counts_folder)

    """
    Create BAM index files.

    """
    print("# experiment BAMs:  %i" % (len(args.exp_bam_list)))
    print("# control BAMs:     %i" % (len(args.ctr_bam_list)))

    exp_bams = []
    control_bams = []
    exp_ids = []
    control_ids = []
    params = "-s -v"
    if args.filter_mode == 2:
        params = "-s -u"
    intermediate_files = []

    print("Create BAM index files ... ")
    for i, bam_file in enumerate(args.exp_bam_list):

        idx = i + 1

        out_bam = args.out_folder + "/%s_exp_rep%i.bam" % (data_id, idx)

        if args.filter_bed:
            print("Filter %s by --filter-bed ... " % (bam_file))
            bam_remove_overlap_region_reads(args.filter_bed, bam_file, out_bam, params=params)
            intermediate_files.append(out_bam)
        else:
            out_bam = bam_file
            # shutil.move(bam_file, out_bam)

        pysam.index(out_bam)
        exp_bams.append(out_bam)
        exp_ids.append("%s_exp_rep%i" % (data_id, idx))

    for i, bam_file in enumerate(args.ctr_bam_list):

        idx = i + 1

        out_bam = args.out_folder + "/%s_ctrl_rep%i.bam" % (data_id, idx)

        if args.filter_bed:
            print("Filter %s by --filter-bed ... " % (bam_file))
            bam_remove_overlap_region_reads(args.filter_bed, bam_file, out_bam, params=params)
            intermediate_files.append(out_bam)
        else:
            out_bam = bam_file
            # shutil.move(bam_file, out_bam)

        pysam.index(out_bam)
        control_bams.append(out_bam)
        control_ids.append("%s_ctrl_rep%i" % (data_id, idx))

    """
    htseq-clip extract crosslink sites.

    htseq-clip extract -i bam/Rbp3_total_rep1.bam -e 1 -s m -o sites/Rbp3_total_rep1.bed

    hce_e : -e [1,2] %i
    hce_s : -s %s
    -s {s,i,d,m,e}, --site {s,i,d,m,e}
                            Crosslink site choices, must be one of: s, i, d, m, e
                            s: start site
                            i: insertion site
                            d: deletion site
                            m: middle site
                            e: end site (default: e).
    hce_g : -g offset nt %i
    hce_q : -q 10
    hce_primary : --primary
    hce_c : -c 4
    hce_m : -m 0
    hce_x : -x 500
    hce_l : -l 10000
    hce_f : -f chr_id file

    """

    print("Extract crosslink sites from BAM files ... ")

    extract_params = " -e %i -s %s -g %i -q %i -c %i -m %i -x %i -l %i" % (args.hce_e, args.hce_s, args.hce_g, args.hce_q, args.hce_c, args.hce_m, args.hce_x, args.hce_l)
    if args.hce_primary:
        extract_params += " --primary "
    if args.hce_f:
        assert os.path.exists(args.hce_f), "--hce-f file \"%s\" not found" % (args.hce_f)
        extract_params += " -f %s " % (args.hce_f)

    # Experiment BAMs.
    exp_beds = []
    for i, bam_file in enumerate(exp_bams):

        idx = i + 1

        out_bed = args.out_folder + "/%s_exp_rep%i.bed" % (data_id, idx)

        check_cmd = "htseq-clip extract -i " + bam_file + " -o " + out_bed + extract_params
        output = subprocess.getoutput(check_cmd)

        print(check_cmd)
        print(output)

        assert os.path.exists(out_bed), "htseq-clip extract -o file \"%s\" not found" % (out_bed)
        exp_beds.append(out_bed)
        intermediate_files.append(out_bed)

    # Control BAMs.
    control_beds = []
    for i, bam_file in enumerate(control_bams):

        idx = i + 1

        out_bed = args.out_folder + "/%s_ctrl_rep%i.bed" % (data_id, idx)

        check_cmd = "htseq-clip extract -i " + bam_file + " -o " + out_bed + extract_params
        output = subprocess.getoutput(check_cmd)

        print(check_cmd)
        print(output)

        assert os.path.exists(out_bed), "htseq-clip extract -o file \"%s\" not found" % (out_bed)
        control_beds.append(out_bed)
        intermediate_files.append(out_bed)

    """
    Count crosslink sites in sliding windows.

    htseq-clip count -i sites/Rbp3_total_rep1.bed -a annotation/Rbp3_uv_total_w75s10.txt -o counts/Rbp3_total_rep1.csv

    """

    print("Count reads overlapping with windows ... ")

    count_params = ""
    if args.hcc_unstranded:
        count_params += " --unstranded"
    # Experiment BAMs.
    for i, bed_file in enumerate(exp_beds):

        idx = i + 1

        out_csv = counts_folder + "/%s_exp_rep%i.csv.gz" % (data_id, idx)
        if args.no_zipper:
            out_csv = counts_folder + "/%s_exp_rep%i.csv" % (data_id, idx)

        check_cmd = "htseq-clip count -i " + bed_file + " -a " + args.in_win_bed + " -o " + out_csv + count_params
        output = subprocess.getoutput(check_cmd)

        print(check_cmd)
        print(output)

        assert os.path.exists(out_csv), "htseq-clip count -o file \"%s\" not found" % (out_csv)
        intermediate_files.append(out_csv)

    # Control BAMs.
    for i, bed_file in enumerate(control_beds):

        idx = i + 1

        out_csv = counts_folder + "/%s_ctrl_rep%i.csv.gz" % (data_id, idx)
        if args.no_zipper:
            out_csv = counts_folder + "/%s_ctrl_rep%i.csv" % (data_id, idx)

        check_cmd = "htseq-clip count -i " + bed_file + " -a " + args.in_win_bed + " -o " + out_csv + count_params
        output = subprocess.getoutput(check_cmd)

        print(check_cmd)
        print(output)

        assert os.path.exists(out_csv), "htseq-clip count -o file \"%s\" not found" % (out_csv)
        intermediate_files.append(out_csv)

    """
    Create an R friendly matrix file.

    htseq-clip createMatrix -i counts/ -b Rbp3 -o counts/Rbp3_uv_total_w75s10_counts.txt.gz

    """

    print("Create R-friendly count matrix file ... ")

    # out_matrix = args.out_folder + "/%s_count_matrix.txt.gz" %(data_id)
    out_matrix = args.out_folder + "/count_matrix.txt.gz"
    if args.no_zipper:
        # out_matrix = args.out_folder + "/%s_count_matrix.txt" %(data_id)
        out_matrix = args.out_folder + "/count_matrix.txt"

    check_cmd = "htseq-clip createMatrix -i " + counts_folder + " -b " + data_id + " -o " + out_matrix
    output = subprocess.getoutput(check_cmd)

    print(check_cmd)
    print(output)

    assert os.path.exists(out_matrix), "htseq-clip createMatrix -o file \"%s\" not found" % (out_matrix)

    """
    Create sample info file for DEWSeq.

    Sample Info file format:
    Sample name	Sample type
    Rbp3_uv_rep1	IP
    Rbp3_uv_rep2	IP
    Rbp3_uv_rep3	IP
    Rbp3_total_rep1	SMI
    Rbp3_total_rep2	SMI
    Rbp3_total_rep3	SMI

    """

    print("Write sample info file ... ")
    sample_info_file = args.out_folder + "/sample_info.txt"
    OUTTAB = open(sample_info_file, "w")
    OUTTAB.write("Sample name\tSample type\n")
    for sid in exp_ids:
        OUTTAB.write("%s\tIP\n" % (sid))
    for sid in control_ids:
        OUTTAB.write("%s\tSMI\n" % (sid))
    OUTTAB.close()

    """
    Delete intermediate files.

    """
    if not args.keep_imf:
        print("Delete intermediate files ... ")
        for imf in intermediate_files:
            if os.path.exists(imf):
                os.remove(imf)

    """
    Report.

    """

    print("")
    print("OUTPUT FILES")
    print("============")

    print("Count matrix file (DEWSeq input file):\n%s" % (out_matrix))
    print("Sample info file (DEWSeq input file):\n%s" % (sample_info_file))
    print("")
