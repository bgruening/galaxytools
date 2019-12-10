#!/usr/bin/env python3


from lib import cliplib
import subprocess
import argparse
import shutil
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


OUTPUT FILES
============
    data_id.model
    data_id.params
if not --disable-cv:
    data_id.cv_results
if not --disable-motifs:
    data_id.sequence_motif
    data_id.sequence_motif.png
    if --str-model:
        data_id.structure_motif
        data_id.structure_motif.png
Temporary:
    data_id.predictions
    data_id.profile


  --opt-set-size int  Hyperparameter optimization set size (taken away from both --pos and --neg) (default: 500)
  --opt-pos str       Positive (= binding site) sequences .fa file for hyperparameter optimization (default: take
                      --opt-set-size from --pos)
  --opt-neg str       Negative sequences .fa file for hyperparameter optimization (default: take --opt-set-size
                      from --neg)
  --min-train int     Minimum amount of training sites demanded (default: 500)
  --disable-cv        Disable cross validation step (default: false)
  --disable-motifs    Disable motif generation step (default: false)
  --gp-output         Print output produced by GraphProt (default: false)
  --str-model         Train a structure model (default: train a sequence model)


EXAMPLE CALLS
=============

python graphprot_train_wrapper.py --pos gp_data/SERBP1_positives.train.fa --neg gp_data/SERBP1_negatives.train.fa --data-id test2 --disable-cv --gp-output --opt-set-size 200 --min-train 400

python graphprot_train_wrapper.py --pos gp_data/SERBP1_positives.train.fa --neg gp_data/SERBP1_negatives.train.fa --data-id test2 --disable-cv --opt-set-size 100 --min-train 200


"""

################################################################################

def setup_argument_parser():
    """Setup argparse parser."""
    help_description = """
    Galaxy wrapper script for GraphProt to train a GraphProt model on 
    a given set of input sequences (positives and negatives .fa). By 
    default a sequence model is trained (due to structure models 
    being much slower to train). Also by default take a portion of 
    the input sequences for hyperparameter optimization (HPO) prior to 
    model training, and run a 10-fold cross validation and motif 
    generation after model training. Thus the following output 
    files are produced: 
    .model model file, .params model parameter file, .png motif files 
    (sequence, or sequence+structure), .cv_results CV results file.
    After model training, predict on positives to get highest whole 
    site and profile scores found in binding sites. Take the median
    score out of these to store in .params file, using it later
    for outputting binding sites or peaks with higher confidence.

    """
    # Define argument parser.
    p = argparse.ArgumentParser(add_help=False,
                                prog="graphprot_train_wrapper.py",
                                description=help_description,
                                formatter_class=argparse.MetavarTypeHelpFormatter)

    # Argument groups.
    p_man = p.add_argument_group("MANDATORY ARGUMENTS")
    p_opt = p.add_argument_group("OPTIONAL ARGUMENTS")

    # Required arguments.
    p_opt.add_argument("-h", "--help",
           action="help",
           help="Print help message")
    p_man.add_argument("--pos",
           dest="in_pos_fa",
           type=str,
           required = True,
           help = "Positive (= binding site) sequences .fa file for model training (option -fasta)")
    p_man.add_argument("--neg",
           dest="in_neg_fa",
           type=str,
           required = True,
           help = "Negative sequences .fa file for model training (option -negfasta)")
    p_man.add_argument("--data-id",
           dest="data_id",
           type=str,
           required = True,
           help = "Data ID (option -prefix)")
    # Additional arguments.
    p_opt.add_argument("--opt-set-size",
           dest="opt_set_size",
           type = int,
           default = 500,
           help = "Hyperparameter optimization set size (taken away from both --pos and --neg) (default: 500)")
    p_opt.add_argument("--opt-pos",
           dest="opt_pos_fa",
           type=str,
           help = "Positive (= binding site) sequences .fa file for hyperparameter optimization (default: take --opt-set-size from --pos)")
    p_opt.add_argument("--opt-neg",
           dest="opt_neg_fa",
           type=str,
           help = "Negative sequences .fa file for hyperparameter optimization (default: take --opt-set-size from --neg)")
    p_opt.add_argument("--min-train",
           dest="min_train",
           type = int,
           default = 500,
           help = "Minimum amount of training sites demanded (default: 500)")
    p_opt.add_argument("--disable-cv",
           dest = "disable_cv",
           default = False,
           action = "store_true",
           help = "Disable cross validation step (default: false)")
    p_opt.add_argument("--disable-motifs",
           dest = "disable_motifs",
           default = False,
           action = "store_true",
           help = "Disable motif generation step (default: false)")
    p_opt.add_argument("--gp-output",
           dest = "gp_output",
           default = False,
           action = "store_true",
           help = "Print output produced by GraphProt (default: false)")
    p_opt.add_argument("--str-model",
           dest = "train_str_model",
           default = False,
           action = "store_true",
           help = "Train a structure model (default: train a sequence model)")
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
    assert cliplib.is_tool("GraphProt.pl"), "GraphProt.pl not in PATH"
    # Check file inputs.
    assert os.path.exists(args.in_pos_fa), "positives .fa file \"%s\" not found" %(args.in_pos_fa)
    assert os.path.exists(args.in_neg_fa), "negatives .fa file \"%s\" not found" %(args.in_neg_fa)
    # Count .fa entries.
    c_pos_fa = cliplib.count_fasta_headers(args.in_pos_fa)
    c_neg_fa = cliplib.count_fasta_headers(args.in_neg_fa)
    assert c_pos_fa, "positives .fa file \"%s\" no headers found" %(args.in_pos_fa)
    assert c_neg_fa, "negatives .fa file \"%s\" no headers found" %(args.in_neg_fa)
    print("# positive .fa sequences:   %i" %(c_pos_fa))
    print("# negative .fa sequences:   %i" %(c_neg_fa))
    # Check additional files.
    if args.opt_pos_fa:
        assert args.opt_neg_fa, "--opt-pos but no --opt-neg given"
    if args.opt_neg_fa:
        assert args.opt_pos_fa, "--opt-neg but no --opt-pos given"
    # If parop .fa files given.
    if args.opt_pos_fa and args.opt_neg_fa:
        c_parop_pos_fa = cliplib.count_fasta_headers(args.opt_pos_fa)
        c_parop_neg_fa = cliplib.count_fasta_headers(args.opt_neg_fa)
        assert c_parop_pos_fa, "--opt-pos .fa file \"%s\" no headers found" %(args.opt_pos_fa)
        assert c_parop_neg_fa, "--opt-neg .fa file \"%s\" no headers found" %(args.opt_neg_fa)
        # Less than 500 for training?? You gotta be kidding.
        assert c_pos_train >= args.min_train, "--pos for training < %i, please provide more (try at least > 1000, the more the better)" %(args.min_train)
        assert c_neg_train >= args.min_train, "--neg for training < %i, please provide more (try at least > 1000, the more the better)" %(args.min_train)
        # Looking closer at ratios.
        pos_neg_ratio = c_parop_pos_fa / c_parop_neg_fa
        if pos_neg_ratio < 0.8 or pos_neg_ratio > 1.25:
            assert 0, "ratio of --opt-pos to --opt-neg < 0.8 or > 1.25 (ratio = %f). Try to keep ratio closer to 1 or better use identical numbers (keep in mind that performance measures such as accuracy or AUROC are not suitable for imbalanced datasets!)" %(pos_neg_ratio)
    else:
        # Define some minimum amount of training sites for the sake of sanity.
        c_pos_train = c_pos_fa - args.opt_set_size
        c_neg_train = c_neg_fa - args.opt_set_size
        # Start complaining.
        assert c_pos_fa >= args.opt_set_size, "# positives < --opt-set-size (%i < %i)" %(c_pos_fa, args.opt_set_size)
        assert c_neg_fa >= args.opt_set_size, "# negatives < --opt-set-size (%i < %i)" %(c_neg_fa, args.opt_set_size)
        assert c_pos_train >= args.opt_set_size, "# positives remaining for training < --opt-set-size (%i < %i)" %(c_pos_train, args.opt_set_size)
        assert c_neg_train >= args.opt_set_size, "# negatives remaining for training < --opt-set-size (%i < %i)" %(c_neg_train, args.opt_set_size)
        # Less than 500?? You gotta be kidding.
        assert c_pos_train >= args.min_train, "# positives remaining for training < %i, please provide more (try at least > 1000, the more the better)" %(args.min_train)
        assert c_neg_train >= args.min_train, "# negatives remaining for training < %i, please provide more (try at least > 1000, the more the better)" %(args.min_train)
        # Looking closer at ratios.
        pos_neg_ratio = c_pos_train / c_neg_train
        if pos_neg_ratio < 0.8 or pos_neg_ratio > 1.25:
            assert 0, "ratio of --pos to --neg < 0.8 or > 1.25 (ratio = %f). Try to keep ratio closer to 1 or better use identical numbers (keep in mind that performance measures such as accuracy or AUROC are not suitable for imbalanced datasets!)" %(pos_neg_ratio)

    """
    Generate parop + train .fa output files for hyperparameter optimization + training.
    
    """
    # Output files for training.
    pos_parop_fa = args.data_id + ".positives.parop.fa"
    neg_parop_fa = args.data_id + ".negatives.parop.fa"
    pos_train_fa = args.data_id + ".positives.train.fa"
    neg_train_fa = args.data_id + ".negatives.train.fa"

    # If parop .fa files given.
    if args.opt_pos_fa and args.opt_neg_fa:
        # Just copy parop and train files.
        cliplib.make_file_copy(args.opt_pos_fa, pos_parop_fa)
        cliplib.make_file_copy(args.opt_neg_fa, neg_parop_fa)
        cliplib.make_file_copy(args.in_pos_fa, pos_train_fa)
        cliplib.make_file_copy(args.in_neg_fa, neg_train_fa)
    else:
        # Generate parop + train .fa files from input .fa files.
        cliplib.split_fasta_into_test_train_files(args.in_pos_fa, pos_parop_fa, pos_train_fa,
                                                  test_size=args.opt_set_size)
        cliplib.split_fasta_into_test_train_files(args.in_neg_fa, neg_parop_fa, neg_train_fa,
                                                  test_size=args.opt_set_size)

    """
    Do the hyperparameter optimization.
    
    """
    print("Starting hyperparameter optimization (-action ls) ... ")
    check_cmd = "GraphProt.pl -action ls -prefix " + args.data_id + " -fasta " + pos_parop_fa + " -negfasta " + neg_parop_fa
    # If sequence model should be trained (default).
    if not args.train_str_model:
        check_cmd += " -onlyseq"
    print(check_cmd)
    output = subprocess.getoutput(check_cmd)
    #assert output, "The following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
    #if args.gp_output:
    #    print(output)
    params_file = args.data_id + ".params"
    assert os.path.exists(params_file), "Hyperparameter optimization output .params file \"%s\" not found" %(params_file)
    # Add model type to params file.
    if args.train_str_model:
        cliplib.echo_add_to_file("model_type: structure", params_file)
    else:
        cliplib.echo_add_to_file("model_type: sequence", params_file)
    # Get parameter string.
    param_string = cliplib.graphprot_get_param_string(params_file)

    """
    Do the model training. (Yowza!)
    
    """
    print("Starting model training (-action train) ... ")
    check_cmd = "GraphProt.pl -action train -prefix " + args.data_id + " -fasta " + pos_train_fa + " -negfasta " + neg_train_fa + " " + param_string
    print(check_cmd)
    output = subprocess.getoutput(check_cmd)
    assert output, "The following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
    if args.gp_output:
        print(output)
    model_file = args.data_id + ".model"
    assert os.path.exists(model_file), "Training output .model file \"%s\" not found" %(model_file)

    """
    Do the 10-fold cross validation.
    
    """
    if not args.disable_cv:
        print("Starting 10-fold cross validation (-action cv) ... ")
        check_cmd = "GraphProt.pl -action cv -prefix " + args.data_id + " -fasta " + pos_train_fa + " -negfasta " + neg_train_fa + " " + param_string + " -model " + model_file
        print(check_cmd)
        output = subprocess.getoutput(check_cmd)
        assert output, "The following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
        if args.gp_output:
            print(output)
        cv_results_file = args.data_id + ".cv_results"
        assert os.path.exists(cv_results_file), "CV output .cv_results file \"%s\" not found" %(cv_results_file)

    """
    Do the motif generation.
    
    """
    if not args.disable_motifs:
        print("Starting motif generation (-action motif) ... ")
        check_cmd = "GraphProt.pl -action motif -prefix " + args.data_id + " -fasta " + pos_train_fa + " -negfasta " + neg_train_fa + " " + param_string + " -model " + model_file
        print(check_cmd)
        output = subprocess.getoutput(check_cmd)
        assert output, "The following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
        if args.gp_output:
            print(output)
        seq_motif_file = args.data_id + ".sequence_motif"
        seq_motif_png_file = args.data_id + ".sequence_motif.png"
        assert os.path.exists(seq_motif_file), "Motif output .sequence_motif file \"%s\" not found" %(seq_motif_file)
        assert os.path.exists(seq_motif_png_file), "Motif output .sequence_motif.png file \"%s\" not found" %(seq_motif_png_file)
        if args.train_str_model:
            str_motif_file = args.data_id + ".structure_motif"
            str_motif_png_file = args.data_id + ".structure_motif.png"
            assert os.path.exists(str_motif_file), "Motif output .structure_motif file \"%s\" not found" %(str_motif_file)
            assert os.path.exists(str_motif_png_file), "Motif output .structure_motif.png file \"%s\" not found" %(str_motif_png_file)

    """
    Do whole site predictions on positive training set.
    
    """
    print("Starting whole site predictions on positive training set (-action predict) ... ")
    check_cmd = "GraphProt.pl -action predict -prefix " + args.data_id + " -fasta " + pos_train_fa + " " + param_string + " -model " + model_file
    print(check_cmd)
    output = subprocess.getoutput(check_cmd)
    assert output, "The following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
    if args.gp_output:
        print(output)
    ws_predictions_file = args.data_id + ".predictions"
    assert os.path.exists(ws_predictions_file), "Whole site prediction output .predictions file \"%s\" not found" %(ws_predictions_file)

    """
    Do profile predictions on positive training set.
    
    """
    print("Starting profile predictions on positive training set (-action predict_profile) ... ")
    check_cmd = "GraphProt.pl -action predict_profile -prefix " + args.data_id + " -fasta " + pos_train_fa + " " + param_string + " -model " + model_file
    print(check_cmd)
    output = subprocess.getoutput(check_cmd)
    assert output, "The following call of GraphProt.pl produced no output:\n%s" %(check_cmd)
    if args.gp_output:
        print(output)
    profile_predictions_file = args.data_id + ".profile"
    assert os.path.exists(profile_predictions_file), "Profile prediction output .profile file \"%s\" not found" %(profile_predictions_file)

    """
    Get 50 % score (median) for .predictions and .profile file.
    For .profile, first extract for each site the maximum score, and then 
    from the list of maximum site scores get the median.
    For whole site .predictions, get the median from the site scores list.
    
    """
    print("Getting .profile and .predictions median scores ... ")

    # Whole site scores median.
    ws_pred_median = cliplib.graphprot_predictions_get_median(ws_predictions_file)
    # Profile top site scores median.
    profile_median = cliplib.graphprot_profile_get_top_scores_median(profile_predictions_file, 
                                                                     profile_type="profile")
    ws_pred_string = "pos_train_ws_pred_median: %f" %(ws_pred_median)
    profile_string = "pos_train_profile_median: %f" %(profile_median)
    cliplib.echo_add_to_file(ws_pred_string, params_file)
    cliplib.echo_add_to_file(profile_string, params_file)
    # Average profile top site scores median for extlr 1 to 10.
    for i in range(10):
        i += 1
        avg_profile_median = cliplib.graphprot_profile_get_top_scores_median(profile_predictions_file,
                                                                             profile_type="avg_profile",
                                                                             avg_profile_extlr=i)
                                                                        
        avg_profile_string = "pos_train_avg_profile_median_%i: %f" %(i, avg_profile_median)
        cliplib.echo_add_to_file(avg_profile_string, params_file)

    print("Script: I'm done.")
    print("Author: Good. Now go back to your file system directory.")
    print("Script: Ok.")


"""
OLD CODE:

Do it for:
1 to 10

    p.add_argument("--ap-extlr",
                   dest="ap_extlr",
                   type = int,
                   default = 5,
                   help = "Define average profile up- and downstream extension for averaging scores to produce the average profile. This is used to get the median average profile score, which will be stored in the .params file to later be used in a prediction setting as a second filter value to get more confident peak regions. NOTE that you have to use the same value in model training and prediction! (default: 5)")


    # Check input parameters.
    if args.disable_opt:
        if args.train_str_model:
            assert args.param_abstraction, "--disable-opt and --str-model is set, but not --abstraction"
        assert args.param_r, "ERROR: --disable-opt is set, but not --R"
        assert args.param_d, "ERROR: --disable-opt is set, but not --D"
        assert args.param_epochs, "ERROR: --disable-opt is set, but not --epochs"
        assert args.param_lambda, "ERROR: --disable-opt is set, but not --lambda"
        assert args.param_bitsize, "ERROR: --disable-opt is set, but not --bitsize"

    p.add_argument("--disable-opt",
                   dest = "disable_opt",
                   default = False,
                   action = "store_true",
                   help = "Disable hyperparameter optimization (HPO) (default: optimize hyperparameters)")
    p.add_argument("--R",
                   dest = "param_r",
                   type = int,
                   default = False,
                   help = "GraphProt model R parameter (default: determined by HPO)")
    p.add_argument("--D",
                   dest = "param_d",
                   type = int,
                   default = False,
                   help = "GraphProt model D parameter (default: determined by HPO)")
    p.add_argument("--epochs",
                   dest = "param_epochs",
                   type = int,
                   default = False,
                   help = "GraphProt model epochs parameter (default: determined by HPO)")
    p.add_argument("--lambda",
                   dest = "param_lambda",
                   type = float,
                   default = False,
                   help = "GraphProt model lambda parameter (default: determined by HPO)")
    p.add_argument("--bitsize",
                   dest = "param_bitsize",
                   type = int,
                   default = False,
                   help = "GraphProt model bitsize parameter (default: determined by HPO)")
    p.add_argument("--abstraction",
                   dest = "param_abstraction",
                   type = int,
                   default = False,
                   help = "GraphProt model RNAshapes abstraction level parameter for training structure models (default: determined by HPO)")
"""



