#!/usr/bin/env python3

import subprocess
import argparse
import os


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

python htsc_create_sliding_windows.py --gff test-data/paper_tus.Synechocystis_pSYSM.gff3 --out test_compare_out --hcw-w 50 --hcw-s 20 --no-zipper

Compare:
diff test_compare_out/windows_mapped_to_ids.txt test-data/windows.exp.txt
diff test_compare_out/windows.bed test-data/windows.exp.bed
diff test_compare_out/annotation.bed test-data/annotation.exp.bed

This corresponds to:
htseq-clip annotation -g test-data/paper_tus.Synechocystis_pSYSM.gff3 -o test-data/annotation.exp.bed
htseq-clip createSlidingWindows -i test-data/annotation.exp.bed -w 50 -s 20 -o test-data/windows.exp.bed
htseq-clip mapToId -a test-data/windows.exp.bed -o test-data/windows.exp.txt

More tests:
python htsc_create_sliding_windows.py --gff test-data/paper_tus.Synechocystis_pSYSM.gff3 --out test_compare_out --hcw-w 50 --hcw-s 20 --no-zipper --hca-unsorted

DEWSeq input files:
test_compare_out/windows_mapped_to_ids.txt

"""


################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Based on genomic annotations GFF file (--gff), create sliding window
    annotations with htseq-clip.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="htsc_create_sliding_windows.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--gff",
                   dest="in_gff",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Annotation file GFF (so far tested with hg38 GENCODE format). Also accepts gff.gz as well")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Results output folder")
    # htseq-clip annotation.
    p.add_argument("--hca-unsorted",
                   dest="hca_unsorted",
                   default=False,
                   action="store_true",
                   help="htseq-clip annotation --unsorted parameter. Use this flag if the GFF file is unsorted (default: False)")
    # htseq-clip createSlidingWindows.
    p.add_argument("--hcw-w",
                   dest="hcw_w",
                   type=int,
                   metavar='int',
                   default=50,
                   help="htseq-clip createSlidingWindows -w parameter. Sliding window size. If unsure, try 75-100 (default: 50)")
    p.add_argument("--hcw-s",
                   dest="hcw_s",
                   type=int,
                   metavar='int',
                   default=20,
                   help="htseq-clip createSlidingWindows -s parameter. Step size for sliding window (default: 20)")
    # More.
    p.add_argument("--no-zipper",
                   dest="no_zipper",
                   default=False,
                   action="store_true",
                   help="Do not gzip output files (default: False)")
    return p


################################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_gff), "--gff file \"%s\" not found" % (args.in_gff)

    # Output folder.
    if not os.path.exists(args.out_folder):
        os.makedirs(args.out_folder)

    """
    1) Flatten annotation.
    htseq-clip annotation -g args.in_gff -o annotation.bed

    -o content example:
    Synechocystis	5	105	TU1@TU1@protein_coding@exon@1/1@TU1:exon0001	2	-
    Synechocystis	576	990	TU2@TU2@protein_coding@exon@1/1@TU2:exon0001	2	+
    Synechocystis	809	909	TU3@TU3@protein_coding@exon@1/1@TU3:exon0001	2	-
    Synechocystis	1531	2150	TU4@TU4@protein_coding@exon@1/1@TU4:exon0001	2	+
    Synechocystis	2150	2701	TU6@TU6@protein_coding@exon@1/1@TU6:exon0001	2	+
    ...

    """

    annot_bed = args.out_folder + "/annotation.bed.gz"
    if args.no_zipper:
        annot_bed = args.out_folder + "/annotation.bed"

    print("Convert --gff to BED ... ")
    check_cmd = "htseq-clip annotation -g " + args.in_gff + " -o " + annot_bed
    if args.hca_unsorted:
        check_cmd += " --unsorted"
    output = subprocess.getoutput(check_cmd)

    print(check_cmd)
    print(output)

    assert os.path.exists(annot_bed), "htseq-clip annotation -o file \"%s\" not found" % (annot_bed)

    """
    2) Create sliding windows.
    htseq-clip createSlidingWindows -i annotation.bed -w hcw_w -s hcw_s -o windows.bed.gz

    -o content example:
    Synechocystis	5	80	TU1@TU1@protein_coding@exon@1/1@TU1:exon0001W00001@1	2	-
    Synechocystis	15	90	TU1@TU1@protein_coding@exon@1/1@TU1:exon0001W00002@2	2	-
    Synechocystis	25	100	TU1@TU1@protein_coding@exon@1/1@TU1:exon0001W00003@3	2	-
    Synechocystis	35	105	TU1@TU1@protein_coding@exon@1/1@TU1:exon0001W00004@4	2	-
    Synechocystis	576	651	TU2@TU2@protein_coding@exon@1/1@TU2:exon0001W00001@1	2	+
    Synechocystis	586	661	TU2@TU2@protein_coding@exon@1/1@TU2:exon0001W00002@2	2	+
    Synechocystis	596	671	TU2@TU2@protein_coding@exon@1/1@TU2:exon0001W00003@3	2	+
    Synechocystis	606	681	TU2@TU2@protein_coding@exon@1/1@TU2:exon0001W00004@4	2	+
    Synechocystis	616	691	TU2@TU2@protein_coding@exon@1/1@TU2:exon0001W00005@5	2	+
    ...

    """

    win_bed = args.out_folder + "/windows.bed.gz"
    if args.no_zipper:
        win_bed = args.out_folder + "/windows.bed"

    print("Create sliding windows BED ... ")
    win_params = " -w %i -s %i " % (args.hcw_w, args.hcw_s)
    check_cmd = "htseq-clip createSlidingWindows -i " + annot_bed + win_params + " -o " + win_bed
    output = subprocess.getoutput(check_cmd)

    print(check_cmd)
    print(output)

    assert os.path.exists(annot_bed), "htseq-clip createSlidingWindows -o file \"%s\" not found" % (win_bed)

    """
    3) Create mapping file for DEWSeq.
    mapToId: extract "name" column from the annotation file and map the entries
    to unique id and print out in tab separated format.
    htseq-clip mapToId -a windows.bed.gz -o windows.txt.gz

    -o content example:
    unique_id	chromosome	begin	end	strand	gene_id	gene_name	gene_type	gene_region	Nr_of_region	Total_nr_of_region	window_number
    TU1:exon0001W00001	Synechocystis	5	80	-	TU1	TU1	protein_coding	exon	1	1	1
    TU1:exon0001W00002	Synechocystis	15	90	-	TU1	TU1	protein_coding	exon	1	1	2
    TU1:exon0001W00003	Synechocystis	25	100	-	TU1	TU1	protein_coding	exon	1	1	3
    TU1:exon0001W00004	Synechocystis	35	105	-	TU1	TU1	protein_coding	exon	1	1	4
    TU2:exon0001W00001	Synechocystis	576	651	+	TU2	TU2	protein_coding	exon	1	1	1
    TU2:exon0001W00002	Synechocystis	586	661	+	TU2	TU2	protein_coding	exon	1	1	2
    TU2:exon0001W00003	Synechocystis	596	671	+	TU2	TU2	protein_coding	exon	1	1	3
    TU2:exon0001W00004	Synechocystis	606	681	+	TU2	TU2	protein_coding	exon	1	1	4
    TU2:exon0001W00005	Synechocystis	616	691	+	TU2	TU2	protein_coding	exon	1	1	5
    ...

    """

    mapped2ids_txt = args.out_folder + "/windows_mapped_to_ids.txt.gz"
    if args.no_zipper:
        mapped2ids_txt = args.out_folder + "/windows_mapped_to_ids.txt"

    print("Create DEWSeq input annotation file ... ")
    win_params = " -w %i -s %i " % (args.hcw_w, args.hcw_s)
    check_cmd = "htseq-clip mapToId -a " + win_bed + " -o " + mapped2ids_txt
    output = subprocess.getoutput(check_cmd)

    print(check_cmd)
    print(output)

    assert os.path.exists(mapped2ids_txt), "htseq-clip mapToId -o file \"%s\" not found" % (mapped2ids_txt)

    """
    Report.

    """

    print("")
    print("OUTPUT FILES")
    print("============")
    print("Annotation BED:\n%s" % (annot_bed))
    print("Windows BED:\n%s" % (win_bed))
    print("Windows mapped to IDs TXT (DEWseq annotation file):\n%s" % (mapped2ids_txt))
    print("")
