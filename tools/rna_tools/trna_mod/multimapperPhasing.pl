#!/usr/bin/env perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2017-08-03 13:12:13 fall>

use strict;
use warnings;
#use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use List::Util 'shuffle';

## define variables
my $sam_fh                    = '';
my $out_fh                    = 'phased.nSorted.sam';
my $mode                      = 'identical';
my $pass_read_end_dist        = 1;
my $pass_read_indel_dist      = 1;
my $verbose                   = 0;
my $choice                    = 'phased';

## define stat variables
my $seen_aln                  = 0;
my $seen_reads                = 0;
my $assigned_reads            = 0;
my $assigned_blocks           = 0;
my $partially_assigned_blocks = 0;
my $omitted_reads             = 0;
my $omitted_blocks            = 0;
my $reads_omitted_qual        = 0;

## get options from Getopt
pod2usage(-verbose => 0)
    unless GetOptions(
		"sam=s"		=> \$sam_fh,
		"out=s"		=> \$out_fh,
		"mode=s"	=> \$mode,
		"ed=i"		=> \$pass_read_end_dist,
		"id=i"		=> \$pass_read_indel_dist,
		"choice=s"  => \$choice,
		"verbose:i" => \$verbose
);


## open filehandler to output
if ( $out_fh ){
  open OUT, "> $out_fh" or die "Error:\tCan't write to file $out_fh\n";
}
else{
  print STDERR "Error:\tno output file name specified; use option -out <FILE>\n";
  die;
}

## read in sam entries line by line and process block of same read
## samfile must be -n sorted
## use <samtools sort -n -O sam test.bam  > test.nSorted.sam> to create it
open SAM, "< $sam_fh" or die "Can't read $sam_fh\n";

my $seen_read_name = 'l71|6Wl_nr650-fXYCf';  # random read name to capture first block occurrence
my @block_bam_lines = ();
my @bam_block_raw   = ();

while(<SAM>){
  
  # skip sam header
  if (m/^@/){
    print OUT $_;
    next;
  };

  #if we want all we print and skip phasing
  if ($choice eq 'all'){
    print OUT $_;
    next;	  
  }
  elsif ($choice eq 'unique'){  #if we want unique we print unique and skip phasing
	  if (m/\tNH:i:1\t/){
		  print OUT $_;
	  }
	  next;
  }
  chomp;

  # count seen alignments
  $seen_aln++;
  
  my @F = split/\t/, $_;
  my $read       = $F[0];
  my $strand     = &parse_strand_bitflag($F[1]);
  my $chromosome = $F[2];
  my $left_start = $F[3];
  my $cigar      = $F[5];
  my $sequence   = $F[9];
  my $right_end  = $left_start + &cigarlength($cigar);
  my $md_flag    = &cigarlength($cigar);

  if ($mode eq "identical"){
	  if(m/MD:Z:([^\s]+)/){
		  $md_flag = $1;
	  }
  }
  elsif($mode eq "noFilter"){
	  $md_flag = "irrelevant";
  }

  if ($read eq "$seen_read_name" || $seen_read_name eq 'l71|6Wl_nr650-fXYCf'){
	  # add to block
	  push @block_bam_lines, [$read, $strand, $chromosome, $left_start, $cigar, $sequence, $right_end, $md_flag];
	  push @bam_block_raw, $_;
	  $seen_read_name =  $read if ($seen_read_name eq 'l71|6Wl_nr650-fXYCf');
  }
  elsif($seen_read_name){
	  # analyse block
	  my $arr_index = &block_analysis(\@block_bam_lines);
	  &block_report(\@bam_block_raw, $arr_index);
	  
	  # delete current block and open new block
	  $seen_read_name =  $read;
	  @block_bam_lines = ();
	  @bam_block_raw   = ();
	  $seen_reads++;
	  push @block_bam_lines, [$read, $strand, $chromosome, $left_start, $cigar, $sequence, $right_end, $md_flag];
	  push @bam_block_raw, $_;
  }
}

# analyse and report last block
my $arr_index = &block_analysis(\@block_bam_lines);
&block_report(\@bam_block_raw, $arr_index);
$seen_reads++;


## print stats
if ($verbose){
  print STDERR "STATS:\tinput alignments:\t$seen_aln\n";
  print STDERR "STATS:\tinput reads:\t\t$seen_reads\n";
  print STDERR "STATS:\tassigned alignments:\t$assigned_reads\n";
  print STDERR "STATS:\tassigned read blocks:\t$assigned_blocks\n";
  print STDERR "STATS:\tpartially assigned read blocks:\t$partially_assigned_blocks\n";
  print STDERR "STATS:\tomitted alignments:\t$omitted_reads\n";
  print STDERR "STATS:\tomitted alignments due to qual:\t$reads_omitted_qual\n";
  print STDERR "STATS:\tomitted read blocks:\t$omitted_blocks\n";
}

###################
##  SUBROUTINES  ##
###################

sub parse_strand_bitflag{
    my $bitflag = shift;
    my $strand  = '.';
    if($bitflag & 16){
	$strand  = '-';
    }
    else{
	$strand  = '+';
    }
    return($strand);
}

sub cigarlength {
  my $cigar_string = shift;
  my $cigar_length = 0;

  while($cigar_string=~m/(\d+)[MDX=N]/g){
    $cigar_length+=$1;
  }
  
  if($cigar_length == 0) {
    print STDERR "WARNING:\tCIGAR string <'$cigar_string'> seems corrupt;\n";
  }
  return($cigar_length);
}

sub block_report{
  my $bam_block_raw = shift;
  my $arr_index     = shift;

  my $block_size = scalar(@$bam_block_raw);
  my $report_block_size = scalar(@$arr_index);

  if($$arr_index[0] == -1){
    $omitted_reads   += $block_size;
    $omitted_blocks  ++;
    if ($verbose == 2){
      foreach my $line (0 .. $block_size-1){
        print STDERR "LOG:\tREMOVED BLOCK ($omitted_blocks)\t$bam_block_raw[$line]\n";
      }
    }
  }
  else{
    $assigned_reads += scalar(@$arr_index);
    $assigned_blocks++;

    if( scalar(@$arr_index) !=  scalar(@bam_block_raw)){
      $partially_assigned_blocks ++;
      $omitted_reads += abs(scalar(@$arr_index) -  scalar(@bam_block_raw));
    }

    foreach my $line (0 .. $block_size-1){
      if (grep {m/^$line$/} @$arr_index){
        print OUT "$bam_block_raw[$line]\n";
      }
      else{
        print STDERR "LOG:\tREMOVED FROM PARTIALLY ACCEPTED BLOCK ($partially_assigned_blocks)\t$bam_block_raw[$line]\n" if ($verbose == 2);
      }
    }
  }
}

sub block_analysis  {
    my $arr_ref = shift;
    my @block   = @$arr_ref;
    
    my %md_flags = ();
    my $read_name       = 'XX';
    
    foreach my $block_line (0 .. $#block) {
	$read_name  = $block[$block_line]->[0];
	my $strand     = $block[$block_line]->[1];
	my $chromosome = $block[$block_line]->[2];
	my $left_pos   = $block[$block_line]->[3];
	my $cigar      = $block[$block_line]->[4];
	my $sequence   = $block[$block_line]->[5];
	my @sequence   = split//, $sequence;
	my $right_pos  = $block[$block_line]->[6];
        my $md_flags   = $block[$block_line]->[7];
        
        my $md_eval    = &md_eval($md_flags);
        if($md_eval == 1){
          $md_flags{$md_flags}->{$block_line}++;
        }
    }
    
    my $chosen_md_flag = "";
    my $chosen_md_flag_count = 0;
    foreach my $flag (keys %md_flags){
      my $count = scalar(keys %{$md_flags{$flag}});
      if ($count > $chosen_md_flag_count){
        $chosen_md_flag = $flag;
        $chosen_md_flag_count = $count;
      }
    }

    if( scalar(keys %md_flags) == 1){
      return ([keys %{$md_flags{$chosen_md_flag}}]);
    }
    else{
      return ([-1]);
    }
}

sub md_eval{
  my $md_flag = shift;
  my $judgment = 1;

  # check if there are no mismatches/indels $pass_read_end_dist from the start of the read
  if($md_flag =~ m/^(\d+)/){
    $judgment = 0 if ($1 < $pass_read_end_dist);
  }

  # check if there are no mismatches/indels $pass_read_end_dist from the end of the read
  if($md_flag =~ m/(\d+)$/){
    $judgment = 0 if ($1 < $pass_read_end_dist);
  }

  # check if there are no mismatches/indels $pass_read_indel_dist around an indel
  if($md_flag =~ m/([0-9]+)\^[A-Z]+([0-9]+)/){
    $judgment = 0 if ($1 < $pass_read_indel_dist || $2 < $pass_read_indel_dist);
  }

  $reads_omitted_qual++ if ($judgment == 0);

  return($judgment);
}

