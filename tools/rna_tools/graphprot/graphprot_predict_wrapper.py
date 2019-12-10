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
    p_man = p.add_argument_group("MANDATORY ARGUMENTS")
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
           choices = [0,1,2,3,4,5,6,7,8,9,10],
           help = "Score threshold for extracting peak regions (default: 0)")
    p_opt.add_argument("--max-merge-dist",
           dest="max_merge_dist",
           type = int,
           default = 0,
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
           choices = [1,2,3,4,5,6,7,8,9,10],
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
    assert cliplib.is_tool("GraphProt.pl"), "GraphProt.pl not in PATH"
    # Check file inputs.
    assert os.path.exists(args.in_fa), "input .fa file \"%s\" not found" %(args.in_fa)
    assert os.path.exists(args.in_model), "input .model file \"%s\" not found" %(args.in_model)
    assert os.path.exists(args.in_params), "input .params file \"%s\" not found" %(args.in_params)
    # Count .fa entries.
    c_in_fa = cliplib.count_fasta_headers(args.in_fa)
    assert c_in_fa, "input .fa file \"%s\" no headers found" %(args.in_fa)
    print("# input .fa sequences:   %i" %(c_in_fa))
    # Read in FASTA sequences to check for uppercase sequences.
    seqs_dic = cliplib.read_fasta_into_dic(args.in_fa)
    c_uc_nt = cliplib.seqs_dic_count_uc_nts(seqs_dic)
    assert c_uc_nt, "no uppercase nucleotides in input .fa sequences. Please change sequences to uppercase (keep in mind GraphProt only scores uppercase regions (according to its viewpoint concept))"
    if not args.ws_pred:
        # Check for lowercase sequences.
        c_lc_nt = cliplib.seqs_dic_count_lc_nts(seqs_dic)
        assert not c_lc_nt, "lowercase nucleotides not allowed in profile predictions, since GraphProt only scores uppercase regions (according to its viewpoint concept))"
    # Check .bed.
    if args.genomic_sites_bed:
        # An array of checks, marvelous.
        assert os.path.exists(args.genomic_sites_bed), "genomic .bed file \"%s\" not found" %(args.genomic_sites_bed)
        # Check .bed for content.
        assert cliplib.count_file_rows(args.genomic_sites_bed), "genomic .bed file \"%s\" is empty" %(args.genomic_sites_bed)
        # Check .bed for 6-column format.
        assert cliplib.bed_check_six_col_format(args.genomic_sites_bed), "genomic .bed file \"%s\" appears to not be in 6-column .bed format" %(args.genomic_sites_bed)
        # Check for unique column 4 IDs.
        assert cliplib.bed_check_unique_ids(args.genomic_sites_bed), "genomic .bed file \"%s\" column 4 IDs not unique" %(args.genomic_sites_bed)
        # Read in .bed regions, compare to FASTA sequences (compare IDs + lengths)
        seq_len_dic = cliplib.get_seq_lengths_from_seqs_dic(seqs_dic)
        reg_len_dic = cliplib.bed_get_region_lengths(args.genomic_sites_bed)
        for seq_id in seq_len_dic:
            seq_l = seq_len_dic[seq_id]
            assert seq_id in reg_len_dic, "sequence ID \"\" missing in input .bed \"\"" %(seq_id, args.genomic_sites_bed)
            reg_l = reg_len_dic[seq_id]
            assert seq_l == reg_l, "sequence length differs from .bed region length (%i != %i)" %(seq_l, reg_l)
    # Read in model parameters.
    param_dic = cliplib.graphprot_get_param_dic(args.in_params)
    # Create GraphProt parameter string.
    param_string = cliplib.graphprot_get_param_string(args.in_params)

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
            cliplib.graphprot_filter_predictions_file(ws_predictions_file, filt_ws_predictions_file,
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
        seq_ids_list = cliplib.fasta_read_in_ids(args.in_fa)
        # Calculate average profiles.
        print("Getting average profile from profile (extlr for smoothing: %i) ... " %(args.ap_extlr))
        avg_prof_file = args.data_id + ".avg_profile"
        cliplib.graphprot_profile_calculate_avg_profile(profile_predictions_file,
                                                        avg_prof_file,
                                                        ap_extlr=args.ap_extlr,
                                                        seq_ids_list=seq_ids_list,
                                                        method=2)
        # Extract peak regions on sequences with threshold score 0.
        avg_prof_peaks_file = args.data_id + ".avg_profile.peaks.bed"
        print("Extracting peak regions from average profile (score threshold = 0) ... ")
        cliplib.graphprot_profile_extract_peak_regions(avg_prof_file, avg_prof_peaks_file,
                                               max_merge_dist=0,
                                               sc_thr=args.score_thr)
        # Convert peaks to genomic coordinates.
        if args.genomic_sites_bed:
            avg_prof_gen_peaks_file = args.data_id + ".avg_profile.genomic_peaks.bed"
            print("Converting peak regions to genomic coordinates ... ")
            cliplib.bed_peaks_to_genomic_peaks(avg_prof_peaks_file, avg_prof_gen_peaks_file,
                                               genomic_sites_bed=args.genomic_sites_bed)
        # Extract peak regions with threshold score p50.
        if args.conf_out:
            sc_id = "pos_train_avg_profile_median_%i" %(args.ap_extlr)
            # Filter by pos_train_ws_pred_median median.
            assert sc_id in param_dic, "average profile extlr %i median information missing in .params file" %(args.ap_extlr)
            p50_sc_thr = float(param_dic[sc_id])
            avg_prof_peaks_p50_file = args.data_id + ".avg_profile.p50.peaks.bed"
            print("Extracting p50 peak regions from average profile (score threshold = %f) ... " %(p50_sc_thr))
            cliplib.graphprot_profile_extract_peak_regions(avg_prof_file, avg_prof_peaks_p50_file,
                                                           max_merge_dist=0,
                                                           sc_thr=p50_sc_thr)
            # Convert peaks to genomic coordinates.
            if args.genomic_sites_bed:
                avg_prof_gen_peaks_p50_file = args.data_id + ".avg_profile.p50.genomic_peaks.bed"
                print("Converting p50 peak regions to genomic coordinates ... ")
                cliplib.bed_peaks_to_genomic_peaks(avg_prof_peaks_p50_file, avg_prof_gen_peaks_p50_file,
                                                   genomic_sites_bed=args.genomic_sites_bed)

    print("Script: I'm done.")
    print("Author: ... ")





"""

    p.add_argument("--seq-uc",
                   dest = "seq_uc",
                   default = False,
                   action = "store_true",
                   help = "Convert input sequences to uppercase (default: false)")

SERBP1_K562_rep01_939	1	1.30147
SERBP1_K562_rep01_2175	1	0.621339
SERBP1_K562_rep01_1370	-1	-0.707584
SERBP1_K562_rep01_544	1	0.975721
SERBP1_K562_rep01_1221	1	1.27005
SERBP1_K562_rep02_996	1	1.4934


    Example .params content:
    epochs: 30
    lambda: 0.001
    R: 1
    D: 4
    bitsize: 14
    model_type: sequence
    pos_train_ws_pred_median: 1.033690
    pos_train_profile_median: 8.680340
    pos_train_avg_profile_median: 4.027981
    avg_profile_extlr: 5






##################################
# RUN GP PROFILE PREDICTION.
##################################

# Read in FASTA file headers.
my @fasta_ids;
my $headers = qx/grep ">" $i_fasta/;
while ($headers =~ />(.+?)\n/g) {
    push(@fasta_ids,$1);
}

# Run GP profile prediction.
my $gp_call = "GraphProt.pl -action predict_profile -model $i_model -fasta $i_fasta -prefix $data_id $params_string";
print STDOUT "GraphProt call: $gp_call\n";
# &> profile_prediction.log
qx/$gp_call/;


####################################
# CALCULATE AVERAGE PROFILE.
####################################

# Calculate .average_profile from GraphProt .profile file.
# Also add p-value column if distr parameters set.
my $add_p = 0;
if ($distr_my and $distr_sigma and $distr_xi) {
    $add_p = 1;
}

# Average_profile: 1-based, FASTA headers as IDs, using ap_extlr for averaging scores.
my $profile = "$data_id.profile";
my $average_profile = "$data_id.average_profile";
# Calculate window size.
my $win = $ap_extlr * 2 + 1;

open (IN, $profile) or die "Cannot open $profile: $!\n";
open (OUT, '>', $average_profile) or die "Cannot open $average_profile: $!";

# Old ID.
my $old_id = "-";
# Current ID.
my $cur_id = "-";
# Start position of the window.
my $pos_inc = 0;
# Score array.
my @scores;
# Input row counter.
my $c_in = 0;
# Output row counter.
my $c_out = 0;

while (<IN>) {

    if ($_ =~ /(.+?)\t\d+?\t(.+?)\n/) {

        $cur_id = $1;
        my $score = $2;
        $c_in++;

        # Case: New refseq ID / new paragraph.
        if ($cur_id ne $old_id) {

            # Print remaining entries at the end of former paragraph.
            while (scalar(@scores) >= ($ap_extlr + 1)) {
                # Calculate avg array score.
                my $mean = array_mean(@scores);
                $c_out++;
                $pos_inc++; 
                if ($add_p) {
                    my $y = exp(-1*(1 + ( ($mean-$distr_my)/$distr_sigma )*$distr_xi)**(-1/$distr_xi));
                    my $p = sprintf("%.8f", (1-$y));
                    print OUT "$fasta_ids[$old_id]\t$pos_inc\t$mean\t$p\n";
                } else {
                    print OUT "$fasta_ids[$old_id]\t$pos_inc\t$mean\n";
                }
                # Remove old score from beginning.
                shift(@scores);
            }

            # Reset for new paragraph.
            $pos_inc = 0;
            @scores = ();
            $old_id = $cur_id;
            # Save first score.
            push(@scores, $score);
            next;

        }

        # Push score in array for average calculation.
        push(@scores, $score);

        # Case: array smaller $win.
        if (scalar(@scores) < $win) {
            # Subcase: array more than half size, start with calculation.
            if (scalar(@scores) >= ($ap_extlr + 1)) {
                # Calculate avg array score.
                my $mean = array_mean(@scores);
                $c_out++;
                $pos_inc++;
                if ($add_p) {
                    my $y = exp(-1*(1 + ( ($mean-$distr_my)/$distr_sigma )*$distr_xi)**(-1/$distr_xi));
                    my $p = sprintf("%.8f", (1-$y));
                    print OUT "$fasta_ids[$cur_id]\t$pos_inc\t$mean\t$p\n";
                } else {
                    print OUT "$fasta_ids[$cur_id]\t$pos_inc\t$mean\n";
                }
            }
            next;
        }

        # Case: array has $win size.
        if (scalar(@scores) == $win) {
            # Calculate avg array score.
            my $mean = array_mean(@scores);
            $c_out++;
            $pos_inc++;
            if ($add_p) {
                my $y = exp(-1*(1 + ( ($mean-$distr_my)/$distr_sigma )*$distr_xi)**(-1/$distr_xi));
                my $p = sprintf("%.8f", (1-$y));
                print OUT "$fasta_ids[$cur_id]\t$pos_inc\t$mean\t$p\n";
            } else {
                print OUT "$fasta_ids[$cur_id]\t$pos_inc\t$mean\n";
            }
            # Remove old score from beginning.
            shift(@scores);
            next;
        }
    }
}

# Print remaining entries at the end of last section.
while (scalar(@scores) >= ($ap_extlr + 1)) {
    # Calculate avg array score.
    my $mean = array_mean(@scores);
    $c_out++;
    $pos_inc++;
    if ($add_p) {
        my $y = exp(-1*(1 + ( ($mean-$distr_my)/$distr_sigma )*$distr_xi)**(-1/$distr_xi));
        my $p = sprintf("%.8f", (1-$y));
        print OUT "$fasta_ids[$cur_id]\t$pos_inc\t$mean\t$p\n";
    } else {
        print OUT "$fasta_ids[$cur_id]\t$pos_inc\t$mean\n";
    }
    # Remove old score from beginning.
    shift(@scores);
}
close IN;
close OUT;

qx/rm -f $profile/;


#########################
# GET PEAK REGIONS.
#########################

my $p50_peaks_bed_out = $data_id . ".peak_regions_p50.bed";
my $peaks_bed_out = $data_id . ".peak_regions.bed";

# If p-values were calculated, use set p-value to get peak regions.
if ($add_p) {
    extract_peak_regions_from_p_values($average_profile, $peaks_bed_out, $max_merge_dist, $thr_p);
    # If p50-out set.
    if ($i_p50_out) {
        # If p50 p-value present, also get peak regions file for this threshold.
        if ($p50_p) {
            extract_peak_regions_from_p_values($average_profile, $p50_peaks_bed_out, $max_merge_dist, $p50_p);
        } else {
            qx/touch $p50_peaks_bed_out/;
        }
    }
} else {
    # If no p-values available, use score threshold for defining peak regions.
    extract_peak_regions_from_scores($average_profile, $peaks_bed_out, $max_merge_dist, $thr_sc);
    # If p50-out set.
    if ($i_p50_out) {
        # If p50 score present, also get peak regions file for this threshold.
        if ($p50_sc) {
            extract_peak_regions_from_scores($average_profile, $p50_peaks_bed_out, $max_merge_dist, $p50_sc);
        } else {
            qx/touch $p50_peaks_bed_out/;
        }
    }
}

exit;



################################################################################
################################################################################
# SUBROUTINES.
################################################################################


sub array_mean {

  my $mean = sum(@_)/@_;
  return sprintf("%.5f", $mean);

}


################################################################################

sub extract_peak_regions_from_p_values {

    my ($average_profile, $peak_regions_bed_file, $max_merge_dist, $max_p) = @_;

    my $old_ref = "no";
    my $in = "N";
    my $ref_id;
    my $start;
    my $end;
    my $best_p = 1;
    my $best_sc = -1000;
    my $region_s;
    my $region_e;
    my $zero_pos;

    my $temp_bed_file = "peak_regions.tmp.bed";

    open (IN, $average_profile) or die "Cannot open $average_profile: $!";
    open (OUT, '>', $temp_bed_file) or die "Cannot open $temp_bed_file: $!";

    while (<IN>) {
        chomp;
        my ($ref, $s, $sc, $p) = (split /\t/)[0,1,2,3];
        # If file has zero-based positions.
        if ($s == 0) {
            $zero_pos = 1;
        }
        # If positions are one-based, make them zero-based.
        unless ($zero_pos) {
            $s -= 1;
        }
        # At transcript ends, if in positive region, write and reset.
        if ($old_ref ne $ref) {
            if ($in eq "Y") {
                print OUT "$ref_id\t$region_s\t$region_e\t$end;$best_sc;$best_p\t0\t+\n";
                $in = "N";
            }
        }
        $old_ref = $ref;
        # Deal with good p-value regions.
        if ($p <= $max_p) {
            # Start of a positive cluster.
            if ($in eq "N") {
                $start = $s;
                $region_s = $s;
                $region_e = $s + 1;
                $end = $s + 1;
                $best_p = $p;
                $best_sc = $sc;
                $ref_id = $ref;
                $in = "Y";
                next;
            # Inside a positive cluster.
            } elsif ($in eq "Y") {
                if ($p < $best_p) {
                    $start = $s;
                    $end = $s + 1;
                    $best_p = $p;
                    $best_sc = $sc;
                    $ref_id = $ref;
                }
                $region_e++;
                next;
            }
        } else {
            # If we were in positive cluster before.
            if ($in eq "Y") {
                print OUT "$ref_id\t$region_s\t$region_e\t$end;$best_sc;$best_p\t0\t+\n";
                $in = "N";
            }
        }
    }
    # After last line processed.
    if ($in eq "Y") {
      print OUT "$ref_id\t$region_s\t$region_e\t$end;$best_sc;$best_p\t0\t+\n";
      $in = "N";
    }
    close IN;
    close OUT;
    # If merge distance zero (i.e. end of one block is -1 from start of next block).
    if ($max_merge_dist == 0) {
        qx/cat $temp_bed_file > $peak_regions_bed_file/;
    } else {
        # Merge nearby regions.
        open(IN, $temp_bed_file) or die "Cannot open $temp_bed_file: $!";
        open(OUT, '>', $peak_regions_bed_file) or die "Cannot open $peak_regions_bed_file: $!";
        # For storing current block stats.
        my $block_chr = 0;
        my ($block_s, $block_e, $block_best_pos, $block_best_p, $block_best_sc);
        while (<IN>) {
            chomp;
            my ($chr, $s, $e, $id) = (split /\t/)[0,1,2,3];
            my ($best_pos, $best_sc, $best_p) = (split /;/, $id);
            if ($chr eq $block_chr) {
                # If $block_e, $s within merge merge.
                if ( ($s - $block_e) <= $max_merge_dist ) {
                    # Update block stats.
                    $block_e = $e;
                    if ($block_best_p > $best_p) {
                        $block_best_p = $best_p;
                        $block_best_sc = $best_sc;
                        $block_best_pos = $best_pos;
                    }
                } else {
                    # If $e outside merge range, print block.
                    print OUT "$block_chr\t$block_s\t$block_e\t$block_best_pos;$block_best_sc;$block_best_p\t0\t+\n";
                    # Store new block.
                    ($block_chr, $block_s, $block_e, $block_best_pos, $block_best_sc, $block_best_p) = ($chr, $s, $e, $best_pos, $best_sc, $best_p);
                }

            } else {
                # If new chromosome, print last block, otherwise it is the first block.
                if ($block_chr) {
                    print OUT "$block_chr\t$block_s\t$block_e\t$block_best_pos;$block_best_sc;$block_best_p\t0\t+\n";
                }
                ($block_chr, $block_s, $block_e, $block_best_pos, $block_best_sc, $block_best_p) = ($chr, $s, $e, $best_pos, $best_sc, $best_p);
            }
        }
        # Print last block.
        if ($block_chr) {
            print OUT "$block_chr\t$block_s\t$block_e\t$block_best_pos;$block_best_sc;$block_best_p\t0\t+\n";
        }
        close OUT;
        close IN;
    }
    qx/rm -f $temp_bed_file/;
}


################################################################################

sub extract_peak_regions_from_scores {

    my ($average_profile, $peak_regions_bed_file, $max_merge_dist, $min_sc) = @_;

    my $old_ref = "no";
    my $in = "N";
    my $ref_id;
    my $start;
    my $end;
    my $best_sc = -1000;
    my $region_s;
    my $region_e;
    my $zero_pos;

    my $temp_bed_file = "peak_regions.tmp.bed";

    open (IN, $average_profile) or die "Cannot open $average_profile: $!";
    open (OUT, '>', $temp_bed_file) or die "Cannot open $temp_bed_file: $!";

    while (<IN>) {
        chomp;
        my ($ref, $s, $sc) = (split /\t/)[0,1,2];
        # If file has zero-based positions.
        if ($s == 0) {
            $zero_pos = 1;
        }
        # If positions are one-based, make them zero-based.
        unless ($zero_pos) {
            $s -= 1;
        }
        # At transcript ends, if in positive region, write and reset.
        if ($old_ref ne $ref) {
            if ($in eq "Y") {
                print OUT "$ref_id\t$region_s\t$region_e\t$end;$best_sc\t0\t+\n";
                $in = "N";
            }
        }
        $old_ref = $ref;
        # Deal with positive regions.
        if ($sc > $min_sc) {
            # Start of a positive cluster.
            if ($in eq "N") {
                $start = $s;
                $region_s = $s;
                $region_e = $s + 1;
                $end = $s + 1;
                $best_sc = $sc;
                $ref_id = $ref;
                $in = "Y";
                next;
            # Inside a positive cluster.
            } elsif ($in eq "Y") {
                if ($sc > $best_sc) {
                    $start = $s;
                    $end = $s + 1;
                    $best_sc = $sc;
                    $ref_id = $ref;
                }
                $region_e++;
                next;
            }
        } else {
            # If we were in positive cluster before.
            if ($in eq "Y") {
                print OUT "$ref_id\t$region_s\t$region_e\t$end;$best_sc\t0\t+\n";
                $in = "N";
            }
        }
    }
    # After last line processed.
    if ($in eq "Y") {
      print OUT "$ref_id\t$region_s\t$region_e\t$end;$best_sc\t0\t+\n";
      $in = "N";
    }
    close IN;
    close OUT;
    # If merge distance zero (i.e. end of one block is -1 from start of next block).
    if ($max_merge_dist == 0) {
        qx/cat $temp_bed_file > $peak_regions_bed_file/;
    } else {
        # Merge nearby regions.
        open(IN, $temp_bed_file) or die "Cannot open $temp_bed_file: $!";
        open(OUT, '>', $peak_regions_bed_file) or die "Cannot open $peak_regions_bed_file: $!";
        # For storing current block stats.
        my $block_chr = 0;
        my ($block_s, $block_e, $block_best_pos, $block_best_sc);
        while (<IN>) {
            chomp;
            my ($chr, $s, $e, $id) = (split /\t/)[0,1,2,3];
            my ($best_pos, $best_sc) = (split /;/, $id);
            if ($chr eq $block_chr) {
                # If $block_e, $s within merge merge.
                if ( ($s - $block_e) <= $max_merge_dist ) {
                    # Update block stats.
                    $block_e = $e;
                    if ($block_best_sc < $best_sc) {
                        $block_best_sc = $best_sc;
                        $block_best_pos = $best_pos;
                    }
                } else {
                    # If $e outside merge range, print block.
                    print OUT "$block_chr\t$block_s\t$block_e\t$block_best_pos;$block_best_sc\t0\t+\n";
                    # Store new block.
                    ($block_chr, $block_s, $block_e, $block_best_pos, $block_best_sc) = ($chr, $s, $e, $best_pos, $best_sc);
                }

            } else {
                # If new chromosome, print last block, otherwise it is the first block.
                if ($block_chr) {
                    print OUT "$block_chr\t$block_s\t$block_e\t$block_best_pos;$block_best_sc\t0\t+\n";
                }
                ($block_chr, $block_s, $block_e, $block_best_pos, $block_best_sc) = ($chr, $s, $e, $best_pos, $best_sc);
            }
        
        }
        # Print last block.
        if ($block_chr) {
            print OUT "$block_chr\t$block_s\t$block_e\t$block_best_pos;$block_best_sc\t0\t+\n";
        }
        close OUT;
        close IN;
    }
    qx/rm -f $temp_bed_file/;
}


################################################################################

"""


