#!/usr/bin/env python3

import subprocess
import argparse
import shutil
import gplib
import gzip
import sys
import os


"""

TOOL DEPENDENCIES
=================

GraphProt 1.1.7
Best install via:
https://anaconda.org/bioconda/graphprot
Tested with: miniconda3, conda 4.7.12


Script: What's my job this time, master?
Author: It'll be a though one.
Script: I take this as a given.
Author: Oh yeah?
Script: ... I'm ready.


OUTPUT FILES
============

    data_id.avg_profile
    data_id.avg_profile.peaks.bed
--conf-out
    data_id.avg_profile.p50.peaks.bed
--gen-site-bed
    data_id.avg_profile.genomic_peaks.bed
--conf-out --gen-site-bed
    data_id.avg_profile.p50.genomic_peaks.bed
--ws-pred
    data_id.predictions
--ws-pred --conf-out
    data_id.predictions
    data_id.p50.predictions


EXAMPLE CALLS
=============

python graphprot_predict_wrapper.py --model test2.model --params test2.params --fasta gp_data/test10_predict.fa --data-id test2pred --gp-output
python graphprot_predict_wrapper.py --model test2.model --params test2.params --fasta gp_data/test10_predict.fa --data-id test2pred --gen-site-bed gp_data/test10_predict.bed
python graphprot_predict_wrapper.py --model test2.model --params test2.params --fasta gp_data/test10_predict.fa --data-id test2pred --gen-site-bed gp_data/test10_predict.bed --conf-out
python graphprot_predict_wrapper.py --model test2.model --params test2.params --fasta gp_data/test10_predict.fa --data-id test2pred --conf-out --ws-pred


python graphprot_predict_wrapper.py --model test-data/test.model --params test-data/test.params --fasta test-data/test_predict.fa --data-id predtest

python graphprot_predict_wrapper.py --model test-data/test.model --params test-data/test.params --fasta test-data/test_predict.fa --data-id predtest --gen-site-bed test-data/test_predict.bed

--gen-site-bed test-data/test_predict.bed

python graphprot_train_wrapper.py --pos test-data/test_positives.train.fa --neg test-data/test_negatives.train.fa --data-id gptest2 --disable-cv --disable-motifs --opt-pos test-data/test_positives.parop.fa --opt-neg test-data/test_negatives.parop.fa


"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Galaxy wrapper script for GraphProt (-action predict and -action 
    predict_profile) to compute whole site or position-wise scores for input 
    FASTA sequences.
    By default, profile predictions are calculated, followed by average 
    profiles computions and peak regions extraction from average profiles.
    If --ws-pred is set, whole site score predictions on input sequences
    will be run instead.
    If --conf-out is set, sites or peak regions with a score >= the median 
    score of positive training sites will be output.
    If --gen-site-bed .bed file is provided, peak regions will be output 
    with genomic coordinates too.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="graphprot_predict_wrapper.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Argument groups.
    p_man = p.add_argument_group("REQUIRED ARGUMENTS")
    p_opt = p.add_argument_group("OPTIONAL ARGUMENTS")

    # Required arguments.
    p_opt.add_argument("-h", "--help",
           action="help",
           help="Print help message")
    p_man.add_argument("--fasta",
           dest="in_fa",
           type=str,
           required = True,
           help = "Sequences .fa file to predict on (option -fasta)")
    p_man.add_argument("--model",
           dest="in_model",
           type=str,
           required = True,
           help = "GraphProt model file to use for predictions (option -model)")
    p_man.add_argument("--params",
           dest="in_params",
           type=str,
           required = True,
           help = "Parameter file for given model")
    p_man.add_argument("--data-id",
           dest="data_id",
           type=str,
           required = True,
           help = "Data ID (option -prefix)")
    # ---> I'm  a conditional argument <---
    p_opt.add_argument("--ws-pred",
           dest = "ws_pred",
           default = False,
           action = "store_true",
           help = "Run a whole site prediction instead of calculating profiles (default: false)")
    # Additional arguments.
    p_opt.add_argument("--sc-thr",
           dest="score_thr",
           type = float,
           default = 0,
           help = "Score threshold for extracting average profile peak regions (default: 0)")
    p_opt.add_argument("--max-merge-dist",
           dest="max_merge_dist",
           type = int,
           default = 0,
           choices = [0,1,2,3,4,5,6,7,8,9,10],
           help = "Maximum merge distance for nearby peak regions (default: report all non-overlapping regions)")
    p_opt.add_argument("--gen-site-bed",
           dest="genomic_sites_bed",
           type=str,
           help = ".bed file specifying the genomic regions of the input .fa sequences. Corrupt .bed information will be punished (default: false)")
    p_opt.add_argument("--conf-out",
           dest="conf_out",
           default = False,
           action = "store_true",
           help = "Output filtered peak regions BED file or predictions file (if --ws-pred) using the median positive training site score for filtering (default: false)")
    p_opt.add_argument("--gp-output",
           dest = "gp_output",
           default = False,
           action = "store_true",
           help = "Print output produced by GraphProt (default: false)")
    p_opt.add_argument("--ap-extlr",
           dest="ap_extlr",
           type = int,
           default = 5,
           choices = [0,1,2,3,4,5,6,7,8,9,10],
           help = "Define average profile up- and downstream extension to produce the average profile. The mean over small sequence windows (window length = --ap-extlr*2 + 1) is used to get position scores, thus the average profile is more smooth than the initial profile output by GraphProt (default: 5)")
    return p


################################################################################

if __name__ == '__main__':

    # Setup argparse.
    parser = setup_argument_parser()
    # Read in command line arguments.
    args = parser.parse_args()
    
    """
    Do all sorts of sanity checking.
    
    """
    # Check for Linux.
    assert "linux" in sys.platform, "please use Linux"
    # Check tool availability.
    assert gplib.is_tool("GraphProt.pl"), "GraphProt.pl not in PATH"
    # Check file inputs.
    assert os.path.exists(args.in_fa), "input .fa file \"%s\" not found" %(args.in_fa)
    assert os.path.exists(args.in_model), "input .model file \"%s\" not found" %(args.in_model)
    assert os.path.exists(args.in_params), "input .params file \"%s\" not found" %(args.in_params)
    # Count .fa entries.
    c_in_fa = gplib.count_fasta_headers(args.in_fa)
    assert c_in_fa, "input .fa file \"%s\" no headers found" %(args.in_fa)
    print("# input .fa sequences:   %i" %(c_in_fa))
    # Read in FASTA sequences to check for uppercase sequences.
    seqs_dic = gplib.read_fasta_into_dic(args.in_fa)
    c_uc_nt = gplib.seqs_dic_count_uc_nts(seqs_dic)
    assert c_uc_nt, "no uppercase nucleotides in input .fa sequences. Please change sequences to uppercase (keep in mind GraphProt only scores uppercase regions (according to its viewpoint concept))"
    if not args.ws_pred:
        # Check for lowercase sequences.
        c_lc_nt = gplib.seqs_dic_count_lc_nts(seqs_dic)
        assert not c_lc_nt, "lowercase nucleotides not allowed in profile predictions, since GraphProt only scores uppercase regions (according to its viewpoint concept))"
    # Check .bed.
    if args.genomic_sites_bed:
        # An array of checks, marvelous.
        assert os.path.exists(args.genomic_sites_bed), "genomic .bed file \"%s\" not found" %(args.genomic_sites_bed)
        # Check .bed for content.
        assert gplib.count_file_rows(args.genomic_sites_bed), "genomic .bed file \"%s\" is empty" %(args.genomic_sites_bed)
        # Check .bed for 6-column format.
        assert gplib.bed_check_six_col_format(args.genomic_sites_bed), "genomic .bed file \"%s\" appears to not be in 6-column .bed format" %(args.genomic_sites_bed)
        # Check for unique column 4 IDs.
        assert gplib.bed_check_unique_ids(args.genomic_sites_bed), "genomic .bed file \"%s\" column 4 IDs not unique" %(args.genomic_sites_bed)
        # Read in .bed regions, compare to FASTA sequences (compare IDs + lengths)
        seq_len_dic = gplib.get_seq_lengths_from_seqs_dic(seqs_dic)
        reg_len_dic = gplib.bed_get_region_lengths(args.genomic_sites_bed)
        for seq_id in seq_len_dic:
            seq_l = seq_len_dic[seq_id]
            assert seq_id in reg_len_dic, "sequence ID \"\" missing in input .bed \"\"" %(seq_id, args.genomic_sites_bed)
            reg_l = reg_len_dic[seq_id]
            assert seq_l == reg_l, "sequence length differs from .bed region length (%i != %i)" %(seq_l, reg_l)
    # Read in model parameters.
    param_dic = gplib.graphprot_get_param_dic(args.in_params)
    # Create GraphProt parameter string.
    param_string = gplib.graphprot_get_param_string(args.in_params)

    """
    Run predictions.
    
    """
    if args.ws_pred:
        # Do whole site prediction.
        print("Starting whole site predictions on input .fa file (-action predict) ... ")
        check_cmd = "GraphProt.pl -action predict -prefix " + args.data_id + " -fasta " + args.in_fa + " " + param_string + " -model " + args.in_model
        output = subprocess.getoutput(check_cmd)
        assert output, "the following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
        if args.gp_output:
            print(output)
        ws_predictions_file = args.data_id + ".predictions"
        assert os.path.exists(ws_predictions_file), "Whole site prediction output .predictions file \"%s\" not found" %(ws_predictions_file)
        if args.conf_out:
            # Filter by pos_train_ws_pred_median median.
            assert "pos_train_ws_pred_median" in param_dic, "whole site top scores median information missing in .params file"
            pos_train_ws_pred_median = float(param_dic["pos_train_ws_pred_median"])
            # Filtered file.
            filt_ws_predictions_file = args.data_id + ".p50.predictions"
            print("Extracting p50 sites from whole site predictions (score threshold = %f) ... " %(pos_train_ws_pred_median))
            gplib.graphprot_filter_predictions_file(ws_predictions_file, filt_ws_predictions_file,
                                                      sc_thr=pos_train_ws_pred_median)
    else:
        # Do profile prediction.
        print("Starting profile predictions on on input .fa file (-action predict_profile) ... ")
        check_cmd = "GraphProt.pl -action predict_profile -prefix " + args.data_id + " -fasta " + args.in_fa + " " + param_string + " -model " + args.in_model
        output = subprocess.getoutput(check_cmd)
        assert output, "the following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
        if args.gp_output:
            print(output)
        profile_predictions_file = args.data_id + ".profile"
        assert os.path.exists(profile_predictions_file), "Profile prediction output .profile file \"%s\" not found" %(profile_predictions_file)
        # Get sequence IDs in order from input .fa file.
        seq_ids_list = gplib.fasta_read_in_ids(args.in_fa)
        # Calculate average profiles.
        print("Getting average profile from profile (extlr for smoothing: %i) ... " %(args.ap_extlr))
        avg_prof_file = args.data_id + ".avg_profile"
        gplib.graphprot_profile_calculate_avg_profile(profile_predictions_file,
                                                        avg_prof_file,
                                                        ap_extlr=args.ap_extlr,
                                                        seq_ids_list=seq_ids_list,
                                                        method=2)
        # Extract peak regions on sequences with threshold score 0.
        avg_prof_peaks_file = args.data_id + ".avg_profile.peaks.bed"
        print("Extracting peak regions from average profile (score threshold = 0) ... ")
        gplib.graphprot_profile_extract_peak_regions(avg_prof_file, avg_prof_peaks_file,
                                               max_merge_dist=0,
                                               sc_thr=args.score_thr)
        # Convert peaks to genomic coordinates.
        if args.genomic_sites_bed:
            avg_prof_gen_peaks_file = args.data_id + ".avg_profile.genomic_peaks.bed"
            print("Converting peak regions to genomic coordinates ... ")
            gplib.bed_peaks_to_genomic_peaks(avg_prof_peaks_file, avg_prof_gen_peaks_file,
                                               genomic_sites_bed=args.genomic_sites_bed)
        # Extract peak regions with threshold score p50.
        if args.conf_out:
            sc_id = "pos_train_avg_profile_median_%i" %(args.ap_extlr)
            # Filter by pos_train_ws_pred_median median.
            assert sc_id in param_dic, "average profile extlr %i median information missing in .params file" %(args.ap_extlr)
            p50_sc_thr = float(param_dic[sc_id])
            avg_prof_peaks_p50_file = args.data_id + ".avg_profile.p50.peaks.bed"
            print("Extracting p50 peak regions from average profile (score threshold = %f) ... " %(p50_sc_thr))
            gplib.graphprot_profile_extract_peak_regions(avg_prof_file, avg_prof_peaks_p50_file,
                                                           max_merge_dist=0,
                                                           sc_thr=p50_sc_thr)
            # Convert peaks to genomic coordinates.
            if args.genomic_sites_bed:
                avg_prof_gen_peaks_p50_file = args.data_id + ".avg_profile.p50.genomic_peaks.bed"
                print("Converting p50 peak regions to genomic coordinates ... ")
                gplib.bed_peaks_to_genomic_peaks(avg_prof_peaks_p50_file, avg_prof_gen_peaks_p50_file,
                                                   genomic_sites_bed=args.genomic_sites_bed)
    # Done.
    print("Script: I'm done.")
    print("Author: ... ")


