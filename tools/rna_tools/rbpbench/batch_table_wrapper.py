#!/usr/bin/env python3

import argparse
import os
import re
import subprocess


###############################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Python wrapper for RBPBench Galaxy wrapper to work with collections of
    input BED files (i.e. to process them with rbpbench batch).
    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="batch_table_wrapper.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Required arguments.
    p.add_argument("-h", "--help",
                   action="help",
                   help="Print help message")
    p.add_argument("--table",
                   dest="in_table",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Input table file with data ID, method ID, RBP ID and file name (Galaxy element identifier in dataset collection) for each to be processed dataset by rbpbench batch")
    p.add_argument("--paths",
                   dest="in_paths",
                   type=str,
                   metavar='str',
                   nargs='+',
                   required=True,
                   help="List of Galaxy BED file paths (--files path1 path2 .. )")
    p.add_argument("--ids",
                   dest="in_ids",
                   type=str,
                   metavar='str',
                   nargs='+',
                   required=True,
                   help="List of Galaxy element identifiers, equal to the BED dataset names in the dataset collection (--ids id1 id2 .. )")
    p.add_argument("--genome",
                   dest="in_genome",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Genomic sequences file (currently supported formats: FASTA)")
    p.add_argument("--out",
                   dest="out_folder",
                   type=str,
                   metavar='str',
                   required=True,
                   help="Batch results output folder")
    # Optional batch arguments.
    p.add_argument("--ext",
                   dest="ext_up_down",
                   type=str,
                   metavar='str',
                   default="0",
                   help="Up- and downstream extension of --in sites in nucleotides (nt). Set e.g. --ext 30 for 30 nt on both sides, or --ext 20,10 for different up- and downstream extension (default: 0)")
    p.add_argument("--motif-db",
                   dest="motif_db",
                   type=int,
                   default=1,
                   choices=[1, 2, 3],
                   help="Motif database to use. 1: human RBP motifs full (259 RBPs, 605 motifs, human_v0.1), 2: human RBP motifs full (low frequencies not rounded, human_v0.1_no_round), 3: human RBP motifs eCLIP (107 RBPs, 316 motifs, human_eclip_v0.1) (default: 1)")
    p.add_argument("--fimo-nt-freqs",
                   dest="fimo_nt_freqs",
                   type=str,
                   metavar='str',
                   default=False,
                   help="Provide FIMO nucleotide frequencies (FIMO option: --bifile) file (default: use internal frequencies file optimized for human transcripts)")
    p.add_argument("--fimo-pval",
                   dest="fimo_pval",
                   type=float,
                   metavar='float',
                   default=0.001,
                   help="FIMO p-value threshold (FIMO option: --thresh) (default: 0.001)")
    p.add_argument("--bed-score-col",
                   dest="bed_score_col",
                   type=int,
                   metavar='int',
                   default=5,
                   help="--in BED score column used for p-value calculations. BED score can be e.g. log2 fold change or -log10 p-value of the region (default: 5)")
    p.add_argument("--unstranded",
                   dest="unstranded",
                   default=False,
                   action="store_true",
                   help="Set if --in BED regions are NOT strand-specific, i.e., to look for motifs on both strands of the provided regions. Note that the two strands of a region will still be counted as one region (change with --unstranded-ct) (default: False)")
    p.add_argument("--unstranded-ct",
                   dest="unstranded_ct",
                   default=False,
                   action="store_true",
                   help="Count each --in region twice for RBP hit statistics when --unstranded is enabled. By default, two strands of one region are counted as one region for RBP hit statistics")
    return p


###############################################################################

if __name__ == '__main__':

    parser = setup_argument_parser()
    args = parser.parse_args()

    assert os.path.exists(args.in_table), "--table file \"%s\" not found" % (args.in_file)
    assert os.path.exists(args.in_genome), "--genome file \"%s\" not found" % (args.in_genome)

    c_paths = len(args.in_paths)
    c_ids = len(args.in_ids)
    assert c_paths == c_ids, "given # paths (--paths) != # ids (--ids) (%i != %i). Please provide one ID for each path" % (c_paths, c_ids)

    """
    Check given paths and IDs.

    """

    # Paths.
    paths_dic = {}
    paths_list = []
    for path in args.in_paths:
        assert os.path.exists(path), "--paths %s file not found" % (path)
        if path not in paths_dic:
            paths_dic[path] = 1
        else:
            assert False, "--paths %s given > 1. Please provide unique paths" % (path)
        paths_list.append(path)

    # IDs
    ids_dic = {}
    ids_list = []
    for id in args.in_ids:
        if id not in ids_dic:
            ids_dic[id] = 1
        else:
            assert False, "--ids \"%s\" given > 1. Please provide unique element identifiers (dataset names) inside the dataset collection, in order to unambiguously assign element ID to file path" % (id)
        ids_list.append(id)

    id2path_dic = {}
    for idx, id in enumerate(ids_list):
        path = paths_list[idx]
        id2path_dic[id] = path

    """
    Read in table.

    Column format:
    rbp_id method_id data_id dataset_name

    """

    comb_ids_dic = {}
    id_collect_dic = {}
    id_collect_dic["rbp_id"] = []
    id_collect_dic["method_id"] = []
    id_collect_dic["data_id"] = []
    id_collect_dic["set_name"] = []
    id_collect_dic["path"] = []  # Galaxy file path.

    print("Read in --table ... ")

    with open(args.in_table) as f:
        for line in f:

            if re.search("^#", line):
                continue

            cols = line.strip().split("\t")

            assert len(cols) == 4, "line in --table with # cols != 4 (%i) encountered:%s" % (len(cols), line)

            rbp_id = cols[0]
            method_id = cols[1]
            data_id = cols[2]
            set_name = cols[3]

            if rbp_id == "rbp_id":
                continue

            comb_id = "%s,%s,%s,%s" % (rbp_id, method_id, data_id, set_name)

            if comb_id not in comb_ids_dic:
                comb_ids_dic[comb_id] = 1
            else:
                assert False, "data combination (\"%s\") appears > 1 in --table file. Please provide unique combinations for rbpbench batch calculation" % (comb_id)

            assert set_name in ids_dic, "given dataset name \"%s\" from --table not part of given --ids. Please provide dataset names present in dataset collection" % (set_name)

            id_collect_dic["rbp_id"].append(rbp_id)
            id_collect_dic["method_id"].append(method_id)
            id_collect_dic["data_id"].append(data_id)
            id_collect_dic["set_name"].append(set_name)
            id_collect_dic["path"].append(id2path_dic[set_name])

    f.closed

    assert id_collect_dic["rbp_id"], "nothing read in from --table. Please provide non-empty table in correct format (columns: rbp_id method_id data_id dataset_name)"

    """
    Construct RBPBench batch call.

    """

    batch_call = "rbpbench batch"
    batch_call += " --out %s" % (args.out_folder)
    batch_call += " --genome %s" % (args.in_genome)
    batch_call += " --ext %s" % (args.ext_up_down)
    batch_call += " --motif-db %i" % (args.motif_db)
    if args.fimo_nt_freqs:
        batch_call += " --fimo-nt-freqs %s" % (args.fimo_nt_freqs)
    batch_call += " --fimo-pval %s" % (str(args.fimo_pval))
    batch_call += " --bed-score-col %i" % (args.bed_score_col)
    if args.unstranded:
        batch_call += " --unstranded"
    if args.unstranded_ct:
        batch_call += " --unstranded-ct"

    rbp_ids = (" ").join(id_collect_dic["rbp_id"])
    method_ids = (" ").join(id_collect_dic["method_id"])
    data_ids = (" ").join(id_collect_dic["data_id"])
    paths = (" ").join(id_collect_dic["path"])

    batch_call += " --rbp-list %s" % (rbp_ids)
    batch_call += " --method-list %s" % (method_ids)
    batch_call += " --data-list %s" % (data_ids)
    batch_call += " --bed %s" % (paths)

    """
    Execute RBPBench batch call.
    """

    print("")
    print("EXECUTING CALL:\n%s" % (batch_call))
    output = subprocess.getoutput(batch_call)
    print("")
    print("RUN OUTPUT:\n%s" % (output))
    print("")
    print("DONE.")
