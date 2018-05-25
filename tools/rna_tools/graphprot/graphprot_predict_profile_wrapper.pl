#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Cwd qw(getcwd abs_path);
use List::Util qw(sum);

=head1 NAME

=head1 SYNOPSIS

Galaxy wrapper script for GraphProt (-action predict_profile) to compute 
the binding profile for a given model on a given set of sequences provided 
in FASTA format. After profile prediction, average profiles get computed, 
scores signified and binding peak regions extracted. The score signification 
is done using the provided fitted GEV parameters, either from .params file 
or manually set. If score threshold is set (-thr-sc), p-value assignment 
will be skipped and set score threshold will be used to extract peak 
regions. NOTE: Additional lines .params file are used to store and get 
GEV parameters, as well as type of model (model_type: sequence|structure).
Also, this wrapper currently works for classification mode only.

PARAMETERS:

    -help|h         display help page
    -fasta          Input FASTA file (option -fasta)
    -model          Input .model file (option -model)
    -params         Input .params file
                    NOTE: uses .params file with additional
                    parameters
                    Manually set parameters (below) will override 
                    found settings in .params file
    -data-id        Data ID (option -prefix)

GraphProt model parameters (by default get from .params file):

    -onlyseq        Set if model is a sequence model
    -R              GraphProt model R parameter
    -D              GraphProt model D parameter
    -epochs         GraphProt model epochs parameter
    -lambda         GraphProt model lambda parameter
    -bitsize        GraphProt model bitsize parameter
    -abstraction    GraphProt model RNAshapes abstraction level 
                    parameter (set for structure models)

Peak region extraction parameters:

    -thr-sc         Score threshold for extracting peak regions
                    By default p-value of 0.05 is used. If no p-value 
                    calculation possible, -thr-sc is used with default: 0
    -thr-p          p-value threshold for extracting peak regions
                    By default, peak regions with p = 0.05 are extracted,
                    as well as p50 score peak regions (if info given)
                    Default: 0.05
    -merge-dist     Maximum merge distance for nearby peak regions
                    Default: report all non-overlapping regions
    -p50-out        Output p50 score filtered peak regions BED file
                    default: false

GEV distribution parameters:

    -distr-my       GEV distribution my parameter for calculating p-values
                    from scores
    -distr-sigma    GEV distrubution sigma parameter for calculating 
                    p-values from scores
    -distr-xi       GEV distribution xi parameter for calculating p-values
                    from scores
    -ap-extlr       Used average profile left right extension for 
                    averaging scores, which were used for distribution 
                    fitting. NOTE: usually a value of 5 was used for 
                    for getting GEV distribution and parameters. If you 
                    choose a different value here, calculated p-values 
                    will be wrong!
                    default : 5


=head1 DISCRIPTION

5) Write manual
6) add output p50 file with NOTE

6) put GP into rna_tools

NOTE:
Additional lines .params file used to store and get gev parameters, as well 
as type of model (model_type: sequence|structure).

Example .params content:
epochs: 20
lambda: 0.001
R: 1
D: 4
bitsize: 14
model_type: sequence
#ADDITIONAL MODEL PARAMETERS
ap_extlr: 5
gev_my: -2.5408
gev_sigma: 1.6444
gev_xi: -0.1383
p50_score: 6.51534 
p50_p_val: 0.0009059744 

=cut

############################
# COMMAND LINE CHECKING.
############################

# Command line argument variables.
my ($i_help, $i_fasta, $i_model, $i_params, $i_data_id, $i_thr_sc, $i_thr_p, $i_max_merge_dist, $i_p50_out, $i_gp_r, $i_gp_d, $i_gp_epochs, $i_gp_lambda, $i_gp_bitsize, $i_gp_abstr, $i_gp_onlyseq, $i_distr_my, $i_distr_sigma, $i_distr_xi, $i_ap_extlr);

# Parse the command line.
GetOptions ( "help|h" => \$i_help,
             "fasta:s" => \$i_fasta,
             "model:s" => \$i_model,
             "params:s" => \$i_params,
             "data-id:s" => \$i_data_id,
             "thr-sc:f" => \$i_thr_sc,
             "thr-p:f" => \$i_thr_p,
             "p50-out" => \$i_p50_out,
             "merge-dist:i" => \$i_max_merge_dist,
             "R:i" => \$i_gp_r,
             "D:i" => \$i_gp_d,
             "epochs:i" => \$i_gp_epochs,
             "lambda:f" => \$i_gp_lambda,
             "bitsize:i" => \$i_gp_bitsize,
             "abstr:i" => \$i_gp_abstr,
             "onlyseq" => \$i_gp_onlyseq,
             "distr-my:f" => \$i_distr_my,
             "distr-sigma:f" => \$i_distr_sigma,
             "distr-xi:f" => \$i_distr_xi,
             "ap-extlr:i" => \$i_ap_extlr ) or pod2usage(1);

# Check.
pod2usage(1) if $i_help;
($i_fasta and $i_model) or pod2usage "ERROR: -fasta, -model are mandatory";
if ($i_distr_my or $i_distr_sigma or $i_distr_xi) {
    ($i_distr_my and $i_distr_sigma and $i_distr_xi) or pod2usage "ERROR: expects all three distribution parameters to be set";
}

#######################
# SET PARAMETERS.
#######################

# Prefix.
my $data_id = "GraphProt";
if ($i_data_id) {
    $data_id = $i_data_id;
}
my $thr_sc = 0;
if ($i_thr_sc) {
    $thr_sc = $i_thr_sc;
}
my $thr_p = 0.05;
if ($i_thr_p) {
    $thr_p = $i_thr_p;
}
my $ap_extlr = 5;
if ($i_ap_extlr) {
    $ap_extlr = $i_ap_extlr;
}
my $max_merge_dist = 0;
if ($i_max_merge_dist) {
    $max_merge_dist = $i_max_merge_dist;
}

# Get parameters from .params file.
my %params;
if ($i_params) {
    open(IN, $i_params) or die "Cannot open $i_params: $!";
    while(<IN>) {
        next if ($_ =~ /^#/);
        if ($_ =~ /(.+):\s(.+)/) {
            $params{$1} = $2;
        }
    }
    close IN;
}

# Create GP parameter string.
my $params_string = "";
# -onlyseq
my $model_type = "structure";
if (exists $params{"model_type"}) {
    if ($params{"model_type"} eq "sequence") {
        $params_string .= " -onlyseq";
        $model_type = "sequence";
    }
} elsif ($i_gp_onlyseq) {
    $params_string .= " -onlyseq";
    $model_type = "sequence";
}
# -R
if ($i_gp_r) {
    $params_string .= " -R $i_gp_r";
} elsif (exists $params{"R"}) {
    my $v = $params{"R"};
    $params_string .= " -R $v";
} else {
    die "ERROR: -R needs to be set";
}
# -D
if ($i_gp_d) {
    $params_string .= " -D $i_gp_d";
} elsif (exists $params{"D"}) {
    my $v = $params{"D"};
    $params_string .= " -D $v";
} else {
    die "ERROR: -D needs to be set";
}
# -epochs
if ($i_gp_epochs) {
    $params_string .= " -epochs $i_gp_epochs";
} elsif (exists $params{"epochs"}) {
    my $v = $params{"epochs"};
    $params_string .= " -epochs $v";
} else {
    die "ERROR: -epochs needs to be set";
}
# -lambda
if ($i_gp_lambda) {
    $params_string .= " -lambda $i_gp_lambda";
} elsif (exists $params{"lambda"}) {
    my $v = $params{"lambda"};
    $params_string .= " -lambda $v";
} else {
    die "ERROR: -lambda needs to be set";
}
# -bitsize
if ($i_gp_bitsize) {
    $params_string .= " -bitsize $i_gp_bitsize";
} elsif (exists $params{"bitsize"}) {
    my $v = $params{"bitsize"};
    $params_string .= " -bitsize $v";
} else {
    die "ERROR: -bitsize needs to be set";
}
# -abstraction
if ($i_gp_abstr) {
    $params_string .= " -bitsize $i_gp_abstr";
    die "ERROR: -abstraction set with -onlyseq" unless ($model_type eq "structure");
} elsif (exists $params{"abstraction"}) {
    my $v = $params{"abstraction"};
    $params_string .= " -abstraction $v";
    die "ERROR: -abstraction set with -onlyseq" unless ($model_type eq "structure");
}

# Distribution parameters.
my ($distr_my, $distr_sigma, $distr_xi);
if ($i_distr_my) {
    $distr_my = $i_distr_my;
} elsif (exists $params{"gev_my"}) {
    $distr_my = $params{"gev_my"};
}
if ($i_distr_sigma) {
    $distr_sigma = $i_distr_sigma;
} elsif (exists $params{"gev_sigma"}) {
    $distr_sigma = $params{"gev_sigma"};
}
if ($i_distr_xi) {
    $distr_xi = $i_distr_xi;
} elsif (exists $params{"gev_xi"}) {
    $distr_xi = $params{"gev_xi"};
}
# Average profile extension parameter.
if (exists $params{"ap_extlr"}) {
    $ap_extlr = $params{"ap_extlr"};
}

print STDOUT "model_type:     $model_type\n";
print STDOUT "ap_extlr:       $ap_extlr\n";

if ($distr_my and $distr_sigma and $distr_xi) {
    print STDOUT "distr_my:       $distr_my\n";
    print STDOUT "distr_sigma:    $distr_sigma\n";
    print STDOUT "distr_xi:       $distr_xi\n";
}
print STDOUT "max_merge_dist: $max_merge_dist\n";

# p50 filter score.
my $p50_sc;
if (exists $params{"p50_score"}) {
    $p50_sc = $params{"p50_score"};
    print STDOUT "p50_score:      $p50_sc\n";
}
# p50 p-value.
my $p50_p;
if (exists $params{"p50_p_val"}) {
    $p50_p = $params{"p50_p_val"};
    print STDOUT "p50_p_val:      $p50_p\n";
}



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



