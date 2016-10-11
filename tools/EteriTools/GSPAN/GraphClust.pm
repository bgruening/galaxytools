package GraphClust;

## GraphClust Library
##
## Author: Steffen Heyne (heyne@informatik.uni-freiburg.de)
##
## This perl module contains some general functions for GraphLClust
##

use strict;
use warnings;
use Cwd qw(abs_path getcwd);
use File::Path;
use FindBin;
use Array::Utils qw(:all);

use lib "$FindBin::Bin";

use Statistics::Descriptive;
use List::Util 'shuffle';
use List::Util qw/ min max /;

##use GraphClust_config;    ## configured external tool paths

use vars qw(%CONFIG);
##*CONFIG = \%GraphClust_config::CONFIG;

require Exporter;

our @ISA    = qw(Exporter);
our @EXPORT = qw(
  SECTION
  SUBSECTION
  readConfigFile
  printConfig
  writeConfig
  read_fasta_file
  read_fasta_with_nonunique_headers
  read_fasta_with_nonunique_headers_meta
  readSubset
  writeSet
  writeSubsetFasta
  system_call
  newick_tree_to_postorder2
  getNodeLeafs
  writeAln
  readAln
  getBasePairs
  mlocarna_center
  evalCenters
  evalCLUSTER
  evalSVECTOR
  printTree
);

our @EXPORT_OK = qw(
  %CONFIG
  $node_sym
);

our $node_sym = "\$\$nodesym";

our $BIN_DIR = abs_path($FindBin::Bin);
################################################################################
## default config

our %CONFIG = (

  ## import config from GraphClust_config.pm
  %CONFIG,

  OPTS_nspdk_centers => "-ensf 5 -oc -usn",
  OPTS_nspdk         => "-R 3 -D 3 -gt DIRECTED",

  OPTS_locarna_paligs => "-p 0.001 --max-diff-am 50 --tau 50 --indel-open -400 --indel -200 --struct-weight 180",
  OPTS_locarna_maligs => "-p 0.001 --max-diff-am 50 --tau 50 --alifold-consensus-dp",
  OPTS_locarna_p_model => "-p 0.001 --max-diff-am 50 --tau 50 --struct-weight 160 --min-bm-prob 0.0005 --min-am-prob 0.0005 --plfold-span 150 --plfold-winsize 200 --temperature 180 --mea-beta 400 --consistency-transformation",
  OPTS_locarna_model => "-p 0.001 --max-diff-am 50 --tau 50 --alifold-consensus-dp",
  OPTS_fasta2shrep_gspan => "-t '3=0,5=80' -M 5 -c 20 -win 40,150 -shift 30 --cue -u --ignore-header --stack --seq-graph-t --seq-graph-alph",
  OPTS_RNAplfold => "--noLP -c 0.0005 -L 150 -W 200 ",
  OPTS_RNAfold   => "--noLP  ",

  GLOBAL_num_clusters => 5,
  GLOBAL_iterations   => 2,
  GLOBAL_hit_blacklist_overlap => 0.2, ## required overlap of cmsearch hit with fragment for blacklisting
  GLOBAL_group_size            => 100,
  GLOBAL_plfold_minlen         => 150,

  nspdk_knn_center   => 20,
  nspdk_nhf          => 400,
  nspdk_nhf_max      => 1000,
  nspdk_nhf_step     => 50,
  nspdk_fcs          => 1,
  nspdk_greylist_num => 0,

  center_subtree_min => 3,
  center_subtree_max => 7,
  center_tree_type   => 3,
  ## 1: get tree via full pairwise locarna aligs,
  ##    final tree based on locarna scores
  ##
  ## 2: get tree via full pairwise locarna aligs,
  ##    transform locarna score matrix into similarity matrix
  ##    and build tree from sim matrix
  ##    2 IS NOT IMPLEMENTED YET!!!
  ##
  ## 3: get tree via nspdk kernel similarity
  ##
  ## for 1 & 2: - "blacklist" scores of overlapping frags

  center_tree_aligs => 1,
  ## 1: use locarna for subtree alignments
  ## 2: use locarnaP for subtree alignments

  center_model_type => 5,
  ## before: do ranking on subtree score, use best subtree as center model
  ## 1: use multiple alignment from subtree as center model
  ## 2: realign subtree with locarna with options 'OPTS_locarna_model'
  ## 3: realign subtree with locarnaP with options 'OPTS_locarna_p_model'
  ## 4: realign subtree with locarnaP and shrink based on relplot signal
  ## 5: cmfinder, try to extend motif with seqs from extended tree pool
  ##    first realign with locarnaP if not already done, use always rel-signal

  center_skip_unstable => 0,

  cm_calibrate    => 0,
  cm_min_bitscore => 20,
  cm_bitscore_sig => 1,
  cm_max_eval     => 0,
  cm_top_only     => 1,

  input_blastclust     => 1,
  input_blastclust_id  => 90,
  input_blastclust_len => 0.9,
  input_seq_min_length => 30,
  input_win_size       => 150,
  input_win_shift      => 50,
  input_add_revcompl   => 0,

  results_merge_cluster_ol => 0.66,
  results_merge_overlap    => 0.51,
  results_top_num          => 15,
  results_min_cluster_size => 5,
  results_partition_type   => "soft",

  evaluate              => 0,
  evaluate_min_overlap  => 0.51,
  evaluate_class0_as_fp => 1,

## old options

# paligs_sig        => 1,
# cm_simple_model   => 0,
# cm_use_rel_signal => 1,
# cm_use_locarnap   => 1,
# center_simple_select => 1,
# center_simple_tree   => 1,
# sigma                => 0.1, ## ok
# svm_name           => "data.MODEL",
# svm_subtree_cutoff => 0.5,
#  evaluate_svm        => 0,
#  blacklist_centers_max => 3,     ## max number of centers to be blacklisted if less than blacklist_max_models were found in last round
#  blacklist_center_frac => 0.6,
#  blacklist_curr_only   => 1,     ## temp blacklist only depend on last round, do not include temp-blacklist from round before
#  blacklist_max_models  => 1,     ## number of models which have to be found in one round that NO blacklist is used for next round
);

################################################################################
## subs

sub SECTION {

  my $name = $_[0];

  print "\n=====================================================================\n";
  print $name. "\n";
  print "=====================================================================\n";
}

sub SUBSECTION {

  my $name = $_[0];

  print "\n---------------------------------------------------------------------\n";
  print $name. "\n";
  print "---------------------------------------------------------------------\n";
}

sub readConfigFile {
  my $conf_file  = $_[0];
  my $force_load = $_[1];

  $force_load = 0 if ( !defined($force_load) );

  open( CONF, "$conf_file" ) or die "Cannot open config file $conf_file! Exit...\n\n";

  while ( my $line = <CONF> ) {
    chomp $line;
    next if ( $line =~ /^#.*/ );
    next if ( $line eq "" );

    my @ent = split( " ", $line );

    if ( @ent > 2 ) {
      my $join_opt = join( " ", @ent[ 1 .. $#ent ] );
      my $opt = $ent[0];
      @ent = ( $opt, $join_opt );
    }

    if ( @ent == 2 ) {
      if ( exists $CONFIG{ $ent[0] } ) {
        $CONFIG{ $ent[0] } = $ent[1];

        #print "set $ent[0]=$ent[1]:\n";
      } elsif ($force_load) {
        $CONFIG{ $ent[0] } = $ent[1];
      } else {
        print "option $ent[0] not recognized!\n";
      }
    } elsif ( @ent == 1 ) {
      if ( exists $CONFIG{ $ent[0] } ) {
        $CONFIG{ $ent[0] } = "";
      }

      #print "option $ent[0] is empty!\n";
    } else {
      print "option $line not recognized!\n";
    }
  }
  return %CONFIG;
}

sub printConfig {
  my $config = $_[0];

  my %conf = %{$config};
  print "\n";
  foreach my $key ( sort keys %conf ) {
    print "OPTION: " . $key . " = " . $conf{$key} . ":\n";
  }
}

sub writeConfig {
  my $file = $_[0];

  open( OUT, ">$file" );

  my $key_len = 30;

  foreach my $key ( grep { $_ =~ /^GLOBAL_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }
  print OUT "\n";

  foreach my $key ( grep { $_ =~ /^OPTS_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }
  print OUT "\n";

  foreach my $key ( grep { $_ =~ /^center_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  foreach my $key ( grep { $_ =~ /^cm_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  foreach my $key ( grep { $_ =~ /^evaluate/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  foreach my $key ( grep { $_ =~ /^input_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  foreach my $key ( grep { $_ =~ /^nspdk_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  foreach my $key ( grep { $_ =~ /^results_/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  print OUT "\n";

  foreach my $key ( grep { $_ =~ /^SGE|^PATH|^VERSION/ } sort keys %CONFIG ) {
    my $num = $key_len - length($key);
    print OUT $key . ( " " x $num ) . $CONFIG{$key} . "\n";
  }

  close(OUT);
}

##################################################################################
# This method parses a standard fasta file and returns the sequences in a hash
# with the name as the key. The name of the sequence is taken from the header,
# which is the first word after the '>' symbol until the first space.
# The first word may contain any symbol (except spaces of course).
# Furthermore, the method deals with multiple lines, and returns a single sequence.
# Input:
#		file		The name of the fasta file
# Output:
#		An array with
#	(1)	A hash reference where the key is the item id and the value is the
#		sequence as a string.
#	(2)	An array reference including the ids in the order they are given in the
#		input file, $file. This information is necessary if you need the exact
#		order, which is not given in the hash.
#	(3) A hash reference of the rest of the header line in the fasta file (after the IDs)
##################################################################################
sub read_fasta_file {
  my ($file,$make_unique) = @_;
  my $FUNCTION = "read_fasta_file in Sequences.pm";

  my $id        = "";
  my $seqstring = "";
  my %fasta     = ();
  my %header    = ();
  my @order     = ();
  my $line      = "";
  my %meta      = ();
  my %seq_meta  = ();
  my $uniq_count = 0;

  open( IN_HANDLE, "<$file" ) || die "ERROR in $FUNCTION: " . "Couldn't open the following file in package Tool," . " sub read_fasta_file: $file\n";


  while ( $line = <IN_HANDLE> ) {
    chomp($line);

    # header (can contain one space after > symbol)
    if ( $line =~ /^\>\s?(\S+)\s*([\S*\s*]*)/ ) {
      if ($id) {
        if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {
#          $uniq_count++;
#          $id .= "_$uniq_count";
#          print "Warning! Make Seq-id unique! now $id\n";
         die "ERROR in $FUNCTION: " . "multiple sequence id '$id', consider using function " . "read_fasta_with_nonunique_headers instead of read_fasta_file";
        }
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        $fasta{$id} = $seqstring;
        $meta{$id}  = {%seq_meta};    ## anonymous hash reference

        $seqstring = "";
        undef(%seq_meta);
      }

      if ($make_unique){
        $uniq_count++;
        $id = $1."_$uniq_count";
      }else{
        $id = $1;
      }

      my $head = $2;
      $head = "" if ( !$head );
      $header{$id} = $head;
      push( @order, $id );
    } elsif ( $line =~ /(.+)\s+(#\S+)\s*$/ && $id ) {

      if ( exists $seq_meta{$2} ) {
        $seq_meta{$2} .= $1;
      } else {
        $seq_meta{$2} = $1;
      }

    } else {
      $seqstring .= $line if ($id);
    }
  }

  if ($id) {
    if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {
      #$uniq_count++;
      #$id .= "_$uniq_count";
      #print "Warning! Make Seq-id unique! now $id\n";
      die "ERROR in $FUNCTION: " . "multiple sequence id '$id', consider using function " . "read_fasta_with_nonunique_headers instead of read_fasta_file";
    }
    $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
    $fasta{$id} = $seqstring;
    $meta{$id}  = +{%seq_meta};    ## anonymous hash reference
    $seqstring  = "";
    undef(%seq_meta);

  }

  return ( \%fasta, \@order, \%header, \%meta );
}

##################################################################################
# This method parses a fasta file and is useful if the header ID is not-unique!!
# It returns the header lines in an array and the sequence lines in the same
# order. It is then up to the user to extract the parts from the header that is
# necessary for the script.
# Furthermore, the method deals with multiple lines, and returns a single sequence
# without line breaks.
# Input:
#		file		The name of the fasta file
# Output:
#	(1)	An array reference for each header line
#	(2) An array reference for each sequence in same order as the headers
##################################################################################
sub read_fasta_with_nonunique_headers_meta {
  my ($file) = @_;
  my $FUNCTION = "read_fasta_file in Sequences.pm";

  my $header    = "";
  my $seqstring = "";
  my @headers   = ();
  my @sequences = ();
  my @meta      = ();
  my %seq_meta  = ();
  open( IN_HANDLE, "<$file" ) || die "ERROR in $FUNCTION:\n" . "Couldn't open the following file in package Tool," . " sub read_fasta_file: $file\n";

  while ( my $line = <IN_HANDLE> ) {
    chomp($line);

    # header (can contain one space after > symbol)
    if ( $line =~ /^\>(.*)/ ) {
      if ($header) {
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        push( @headers,   $header );
        push( @sequences, $seqstring );
        push( @meta,      +{%seq_meta} );    ## anonymous hash reference
                                             #print keys %seq_meta;
        $seqstring = "";
        undef(%seq_meta);
      }
      $header = $1;

    } elsif ( $line =~ /(.+)\s+(#\S+)$/ && $header ) {

      if ( exists $seq_meta{$2} ) {
        $seq_meta{$2} .= $1;
      } else {
        $seq_meta{$2} = $1;
      }

    } elsif ($header) {
      $seqstring .= $line
    }
  }

  if ($header) {
    $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
    push( @headers,   $header );
    push( @sequences, $seqstring );
    push( @meta,      +{%seq_meta} );
    $seqstring = "";
    %seq_meta  = ();
  }
  return ( \@headers, \@sequences, \@meta );
}

##################################################################################
# This method parses a fasta file and is useful if the header ID is not-unique!!
# It returns the header lines in an array and the sequence lines in the same
# order. It is then up to the user to extract the parts from the header that is
# necessary for the script.
# Furthermore, the method deals with multiple lines, and returns a single sequence
# without line breaks.
# Input:
#		file		The name of the fasta file
# Output:
#	(1)	An array reference for each header line
#	(2) An array reference for each sequence in same order as the headers
##################################################################################
sub read_fasta_with_nonunique_headers {
  my ($file) = @_;
  my $FUNCTION = "read_fasta_file in Sequences.pm";

  my $header    = "";
  my $seqstring = "";
  my @headers   = ();
  my @sequences = ();
  open( IN_HANDLE, "<$file" ) || die "ERROR in $FUNCTION:\n" . "Couldn't open the following file in package Tool," . " sub read_fasta_file: $file\n";

  while ( my $line = <IN_HANDLE> ) {
    chomp($line);

    # header (can contain one space after > symbol)
    if ( $line =~ /^\>(.*)/ ) {
      if ($header) {
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        push( @headers,   $header );
        push( @sequences, $seqstring );
        $seqstring = "";
      }
      $header = $1;
    } else {
      $seqstring .= $line if ($header);
    }
  }

  if ($header) {
    $seqstring =~ s/\s*//g;        ## do not allow spaces in sequence
    push( @headers,   $header );
    push( @sequences, $seqstring );
    $seqstring = "";
  }
  return ( \@headers, \@sequences );
}

## read line $sub_idx from file, assume line contains ids <space>-seperated, id-set is returned
sub readSubset {
  my $file     = $_[0];
  my $sub_idx  = $_[1];
  my $max_rank = $_[2];

  open( IN, "$file" ) or die "Cannot open file $file! Exit...\n\n";

  my @set = ();

  my $idx = 0;
  while ( my $line = <IN> ) {
    $idx++;
    next if ( $idx != $sub_idx );
    chomp $line;
    @set = split( " ", $line );
    last;
  }
  close(IN);

  print "(GraphClust::readSubset) File $file contains less than $sub_idx lines! \n"
    if ( !@set );

  if (@set) {
    $max_rank = -1 if ( !$max_rank );
    $max_rank = @set if ( $max_rank == -1 || $max_rank > @set );

    @set = @set[ 0 .. ( $max_rank - 1 ) ];
  }
  return ( \@set );
}

## read first $max_rank seqs from $fasta_ref, order is given in $ids_aref
## if max_rank = -1 then take all ids from ids_aref
#sub getSubset {
#  my $ids_aref  = $_[0];
#  my $fasta_ref = $_[1];
#  my $max_rank  = $_[2];
#
#  my %subset = ();
#
#  $max_rank = @{$ids_aref} if ( $max_rank == -1 );
#
#  ## each $id is a array ref with (idx,start,end)
#  foreach my $idx ( 0 .. ( $max_rank - 1 ) ) {
#    last if ( $idx >= @{$ids_aref} );
#
#    my $id = $ids_aref->[$idx];
#    die "Error! $id does not exists in fasta file! Exit...\n\n" if ( !exists $fasta_ref->[0]->{$id} );
#    $subset{$id} = $fasta_ref->[0]->{$id};
#  }
#
#  return \%subset;
#}

sub writeSet {
  my $set_aref = $_[0];
  my $outfile  = $_[1];

  open( OUT, ">$outfile" );
  ## old with sorting
  ##print OUT join( " ", sort { $a <=> $b } @{$set_aref} ) . "\n";
  ## new no sorting
  print OUT join( " ", @{$set_aref} ) . "\n";
  close(OUT);
}

sub writeSubsetFrags {
  my $fa_aref   = $_[0];
  my $set_aref  = $_[1];
  my $outfile   = $_[2];
  my $seqprefix = $_[3];

  open( OUT, ">$outfile" ) or die "Cannot open file $outfile! Exit...\n\n";

  foreach my $key ( sort { $a <=> $b } @{$set_aref} ) {

    die "Error! Seq-ID $key does not esist in fasta! Exit...\n\n" if ( !exists $fa_aref->[0]->{$key} );
    die "Error! Seq-ID $key is not at pos $key in fasta! Exit...\n\n" if ( $key != $fa_aref->[1]->[ $key - 1 ] );
    my $head = $fa_aref->[2]->{$key};
    my $frag = "MISS";
    $frag = $1 if ( $head =~ /($seqprefix\d+#\d+#\d+#.)\s*.*$/ );
    die "FASTA header wrong! Fragment key ($seqprefix) (eg. 'SEQ123#45#67#+')  not found in header $head! Exit...\n\n"
      if ( $frag eq "MISS" );
    print OUT "$key " . $frag . "\n";
  }
  close(OUT);
}

sub writeSubsetFasta {
  my $fa_aref  = $_[0];
  my $ids_aref = $_[1];
  my $out_file = $_[2];
  my $info     = $_[3];

  #  my $seqs_href   = $_[0];
  #  my $out_file    = $_[1];
  #  my $header_href = $_[2];

  open( OUT, ">$out_file" ) or die "Cannot open file $out_file! Exit...\n\n";

  ## order is important !!!
  foreach my $key ( sort { $a <=> $b } @{$ids_aref} ) {

    die "Error! Seq-ID $key does not esist in fasta! Exit...\n\n" if ( !exists $fa_aref->[0]->{$key} );

    print OUT ">$key " . $fa_aref->[2]->{$key} . "\n";
    print OUT $fa_aref->[0]->{$key} . "\n";

    next if ( !$info );
    ## write meta info: given structure
    print OUT $fa_aref->[3]->{$key}->{"#FS"} . " #FS\n" if ( exists $fa_aref->[3]->{$key}->{"#FS"} );
    print OUT $fa_aref->[3]->{$key}->{"#S"} . " #S\n" if ( exists $fa_aref->[3]->{$key}->{"#S"} );
  }

  ## order is important !!!
  #  foreach my $key ( sort { $a <=> $b } keys %{$seqs_href} ) {
  #    print OUT ">$key";
  #    print OUT " " . $header_href->{$key} if ($header_href);
  #    print OUT "\n";
  #
  #    print OUT $seqs_href->{$key} . "\n";
  #  }

  close(OUT);
  return $out_file;
}

## function that extracts a submatrix from $mat_file
## submatrix is given with array-ref $idx, gives rows, idx is 0 based
## than calc an row avg score (sum-all-rows/rows)
sub getMatrixSum {
  my $mat_file = $_[0];
  my $idx      = $_[1];

  open( IN, "$mat_file" );
  my @rows = <IN>;
  @rows = map { my @a = split( " ", $rows[$_] ); $rows[$_] = \@a } 0 .. $#rows;
  close(IN);

  @{$idx} = sort { $a <=> $b } @{$idx};

  my @row_sel = ();
  my $mode    = 1;
  if ( $mode == 1 ) {
    ## extract submatrix for all idx
    @row_sel = map { my $a = $rows[$_]; my @b = @{$a}[ @{$idx} ]; \@b } @{$idx};
  } else {
    ## extract full lines for all idx
    @row_sel = map { my $a = $rows[$_]; $a } @{$idx};
  }

  #print "sel idx:" . join( ":", @{$idx} ), "\n";
  my $sum_all = 0;
  my $num     = @{$idx};
  foreach my $id ( 0 .. $#{$idx} ) {

    print "row idx" . $id . " val:" . join( ":", @{ $row_sel[$id] } ) . " ";

    my $sum = 0;
    map { $sum += $_ } @{ $row_sel[$id] };

    print " sum $sum\n";
    $sum_all += $sum;
  }
  my $avg = $sum_all / ($num);

  #print "matrix avg $avg\n";
  return $avg;
}

sub system_call {
  my $sysCommand = $_[0];
  my $verbose    = $_[1];
  my $time_file  = $_[2];

  if ($verbose) {
    print "--------------------------------------------------------------------------------\n";
    print "CALL $sysCommand\n";
    print "--------------------------------------------------------------------------------\n";
  }

  if ( !$time_file ) {
    system("$sysCommand") == 0
      or die "ERROR ($0) Could not call system command:\n\n$sysCommand\n\n";
  } else {

    ## check if we have gnu time command, otherwise dont use it (e.g. on MAC OSX)
    my $time_failed = 0;
    system("\\time --version 1>/dev/null 2>/dev/null") == 0 or $time_failed = 1;

    if (!$time_failed){
        system("\\time -p -o $time_file --format \"total 1 %e %e user 1 %U %U\" bash -c '$sysCommand' ") == 0
          or die "ERROR ($0) Could not call system command:\n\n$sysCommand\n\n";
    } else{
      system("echo \"total 1 0.0 0.0 user 1 0.0 0.0\" > $time_file ");
      system("bash -c '$sysCommand' ") == 0
        or die "ERROR ($0) Could not call system command:\n\n$sysCommand\n\n";
    }
  }
  return 1;
}

########################################ä
## newick_tree_to_postorder($tree string)
##
## Translates a newick tree format string into a list of nodes in postorder
##
## $tree string in newick tree format
##
## returns list of nodes/leaves in postorder, use $node_sym for inner nodes
##
## dies if tree is not parsable
##
########################################
sub newick_tree_to_postorder2 {
  my ($tree) = @_;

  my @list;

  #$tree =~ s/:[\d.e-]+//g;

  $tree =~ s/;$//;    # allow that the tree string is terminated by ';'

  my $brcount = 0;

  my $preID = -1;
  my @last  = ();
  my %tree  = ();

  for ( my $i = 0 ; $i < length $tree ; $i++ ) {
    my $c = substr $tree, $i, 1;

    if ( $c eq "(" ) {

      $brcount++;
      $preID++;
      push( @last, $preID );

      # print "open preID $preID\n";

      $tree{$preID}->{LEAF}   = 0;
      $tree{$preID}->{NLEAFS} = 0;
      $tree{$preID}->{CHILDS} = [];
      $tree{$preID}->{NAMES}  = [];

    } elsif ( $c eq ")" ) {

      my $l = pop(@last);

      if ($l) {

        #print "close preID $preID last $l\n";
        push( @{ $tree{ $last[$#last] }->{CHILDS} }, $l );
        push( @{ $tree{ $last[$#last] }->{NAMES} },  @{ $tree{$l}->{NAMES} } );
        $tree{ $last[$#last] }->{NLEAFS} += @{ $tree{$l}->{NAMES} };
      }

      ## parse branch length of inner node
      my $rest = substr $tree, $i + 1;

      my $len = 0;
      my $t   = "";
      if ( $rest =~ /^([^(),]+)/ ) {
        $t = $1;
        $t =~ /\:(.*)/;
        $len = $1;
      }
      $tree{$l}->{BRANCH} = $len;
      $i += ( length $t );

      $brcount--;

      if ( $brcount < 0 ) {
        die "Parse error in tree.";
      }

      push @list, $node_sym;
    } elsif ( $c eq "," ) {
      ## ignore, although we could do syntax checking
    } else {
      my $rest = substr $tree, $i;
      $rest =~ /^([^(),]+)/;
      my $t = $1;

      my @part  = split( ":", $t );
      my $token = $part[0];
      my $len   = $part[1];

      if ( !$token ) {
        $t =~ /\:(.*)/;
        $len = $1;
      }

      #print "t $t $len last $last[$#last]\n";

      if ($token) {

        $preID++;

        $tree{$preID}->{LEAF}   = 1;
        $tree{$preID}->{NAME}   = $token;
        $tree{$preID}->{BRANCH} = $len;

        #print "token $token preID $preID last " . $last[$#last] . "\n";
        push( @{ $tree{ $last[$#last] }->{CHILDS} }, $preID );
        push( @{ $tree{ $last[$#last] }->{NAMES} },  $token );
        $tree{ $last[$#last] }->{NLEAFS}++;
      }

      #print "part $t token<".$token."> len<$len>\n";

      $i += ( length $t ) - 1;
      push @list, $token if ($token);
    }
  }

  # printTree( \%tree );
  return ( \@list, \%tree );
}

sub printTree {
  my $tree = $_[0];

  foreach my $node ( sort { $a <=> $b } keys %{$tree} ) {

    if ( $tree->{$node}->{LEAF} ) {
      print "LEAF $node name " . $tree->{$node}->{NAME} . " len " . $tree->{$node}->{BRANCH} . "\n";
    } else {
      print "NODE $node nleafs "
        . $tree->{$node}->{NLEAFS}
        . " childs "
        . join( ":", @{ $tree->{$node}->{CHILDS} } )
        . " names "
        . join( ":", @{ $tree->{$node}->{NAMES} } ) . " len "
        . $tree->{$node}->{BRANCH} . "\n";
    }
  }
}

sub getNodeLeafs {
  my $tree   = $_[0];
  my $nodeID = $_[1];

  my @leafs;    ## leafs in postorder
  my @nodes;    ## postorder of all internal/leaf nodes

  ## inner nodes
  if ( !$tree->{$nodeID}->{LEAF} ) {
    my @childs = @{ $tree->{$nodeID}->{CHILDS} };

    my ( $childLeafs, $childNodes ) = getNodeLeafs( $tree, $childs[0] );
    push( @leafs, @{$childLeafs} );
    push( @nodes, @{$childNodes} );

    ( $childLeafs, $childNodes ) = getNodeLeafs( $tree, $childs[1] );
    push( @leafs, @{$childLeafs} );
    push( @nodes, @{$childNodes} );
    push( @nodes, $nodeID );          # if (!$tree->{$child}->{LEAF});
  }

  ## leafs
  if ( $tree->{$nodeID}->{LEAF} ) {
    return ( [ $tree->{$nodeID}->{NAME} ], [$nodeID] );
  }

  return ( \@leafs, \@nodes );
}

########################################ä
## tree_partitions($tree_postorder ref of list)
##
## generate partitions out of the postorder tree
##
## $tree_postorder ref of list representing tree in postorder (as
## generated by newick_tree_to_postorder)
##
## returns ref of list of refs of partitions of leaves due to the tree; a
## partition is represented as a list/subset of leaves
##
########################################
#sub tree_partitions($) {
#  my $tree_postorder = shift;
#
#  my @result;
#
#  my @stack;
#
#  for my $item (@$tree_postorder) {
#    if ( $item eq $node_sym ) {
#      my @op1 = @{ $stack[-2] };
#      my @op2 = @{ $stack[-1] };
#
#      $#stack -= 2;
#
#      my @op12 = ( @op1, @op2 );
#
#      push @stack, \@op12;
#    } else {
#      push @stack, [$item];
#    }
#    push @result, $stack[-1];
#  }
#  $#result -= 2;    # the last is empty, the one before symmetric
#  return \@result;
#}

## reads aln and writes it bmreliabilities file to generate rel-profile later

sub evalCenters {
  my $center_file = $_[0];
  my $lens_file   = $_[1];
  my $class_file  = $_[2];
  my $out_name    = $_[3];
  my $center_knn  = $_[4];

  #my $class_occ   = $_[5];

  my $r_path = $CONFIG{PATH_R};

  open( LEN, $lens_file );
  my @lens = <LEN>;
  close(LEN);
  chomp(@lens);

  my $eval_class = 1;
  open( CLASS, $class_file ) or $eval_class = 0;
  my @class = <CLASS> if ($eval_class);
  close(CLASS);
  chomp(@class);

  my %classes = ();
  foreach my $idx ( 0 .. $#class ) {

    $class[$idx] =~ /(\d+)\s(\d+)/;

    #print $class[$idx]."\t".$2."\n";
    $class[$idx] = $2;

    $classes{ $class[$idx] }                   = 0;
    $classes{ $class[$idx] . "#KNN" }          = 0;
    $classes{ $class[$idx] . "#CENT" }         = 0;
    $classes{ $class[$idx] . "#CENT_OCC" }     = 0;
    $classes{ $class[$idx] . "#CENT_UNQ" }     = 0;
    $classes{ $class[$idx] . "#CENT_OCC_UNQ" } = 0;
    $classes{ $class[$idx] . "#CENT_OCC_MAX" } = 0;

  }

  my $idx = 0;
  open( CENTERS, $center_file );

  my %centers = ();
  my @set;
  while ( my $line = <CENTERS> ) {
    $idx++;
    chomp $line;
    @set = split( " ", $line );
    $centers{$idx}->{IDS} = \@set;

    #$centers{$idx}->{IDS_KNN} = \@set[ 0 .. ( $center_knn - 1 ) ];

    my @lens_c  = ();
    my @class_c = ();
    foreach my $id (@set) {
      push( @lens_c,  $lens[ $id - 1 ] );
      push( @class_c, $class[ $id - 1 ] );
    }
    $centers{$idx}->{LENS}  = \@lens_c;
    $centers{$idx}->{CLASS} = \@class_c;

  }
  close(CENTERS);
  $center_file =~ /.*\/([^\/]+)$/;
  my $center_file_name = $1;
  my $num_centers      = $idx;

  ## count class occ for centerfile, each centerfile line has similar class occ
  foreach my $id (@set) {
    $classes{ $class[ $id - 1 ] }++ if ($eval_class);
  }

  open( OUT, ">$out_name.lens" );
  foreach my $key ( 1 .. ( keys %centers ) ) {
    print OUT join( " ", @{ $centers{$key}->{LENS} }[ 0 .. ( $center_knn - 1 ) ] ) . "\n";
  }
  close(OUT);

  open( OUT, ">$out_name.lens_all" );
  foreach my $key ( 1 .. ( keys %centers ) ) {
    print OUT join( " ", @{ $centers{$key}->{LENS} } ) . "\n";
  }
  close(OUT);

  if ($eval_class) {
    foreach my $key ( 1 .. ( keys %centers ) ) {
      my %class_occ_cent = ();
      foreach my $cl ( @{ $centers{$key}->{CLASS} }[ 0 .. ( $center_knn - 1 ) ] ) {
        $classes{ $cl . "#KNN" }++;
        $class_occ_cent{$cl} = 0 if ( !exists $class_occ_cent{$cl} );
        $class_occ_cent{$cl}++;
      }

      my @cent_sorted = sort { $class_occ_cent{$b} <=> $class_occ_cent{$a} } keys %class_occ_cent;
      my @cent_sorted_unq = sort { $class_occ_cent{$b} <=> $class_occ_cent{$a} } grep { $_ != 0 } keys %class_occ_cent;

#    print "check cent:$key cl:" . $cent_sorted_unq[ 0 ] . " idx:" . ( 0 ) . " occ:" . $class_occ_cent{ $cent_sorted_unq[ 0 ] } . " 1:".$class_occ_cent{ $cent_sorted_unq[ 1 ] } ." cl1:".$cent_sorted_unq[ 1 ]."\n";
      if ( $class_occ_cent{ $cent_sorted_unq[0] } != $class_occ_cent{ $cent_sorted_unq[1] } ) {

        $classes{ $cent_sorted_unq[0] . "#CENT_UNQ" }++;
        $classes{ $cent_sorted_unq[0] . "#CENT_OCC_UNQ" } += $class_occ_cent{ $cent_sorted_unq[0] };

        print "uniq cent:$key cl:" . $cent_sorted_unq[0] . " idx:" . (0) . " occ:" . $class_occ_cent{ $cent_sorted_unq[0] } . " 1:" . $class_occ_cent{ $cent_sorted_unq[1] } . "\n";
      }

      my $idx = 1;
      while ( ( $idx == 1 ) || ( $idx == 2 && $cent_sorted[0] == 0 ) || ( $idx > 2 && $class_occ_cent{ $cent_sorted[ $idx - 2 ] } == $class_occ_cent{ $cent_sorted[ $idx - 1 ] } ) ) {
        $classes{ $cent_sorted[ $idx - 1 ] . "#CENT" }++;
        $classes{ $cent_sorted[ $idx - 1 ] . "#CENT_OCC" } += $class_occ_cent{ $cent_sorted[ $idx - 1 ] };

        if ( $class_occ_cent{ $cent_sorted[ $idx - 1 ] } > $classes{ $cent_sorted[ $idx - 1 ] . "#CENT_OCC_MAX" } ) {
          $classes{ $cent_sorted[ $idx - 1 ] . "#CENT_OCC_MAX" } = $class_occ_cent{ $cent_sorted[ $idx - 1 ] };
        }

        print "cent:$key cl:" . $cent_sorted[ $idx - 1 ] . " idx:" . ( $idx - 1 ) . " occ:" . $class_occ_cent{ $cent_sorted[ $idx - 1 ] } . "\n";
        $idx++;
      }
    }

    open( ALL, ">$out_name.classes" );
    foreach my $key ( sort { $a <=> $b } grep { !/\#/ } keys %classes ) {
      next if ( $key =~ /#KNN/ || $key =~ /#CENT/ );

      my $avg = 0;
      $avg = $classes{ $key . "#CENT_OCC" } / $classes{ $key . "#CENT" } if ( $classes{ $key . "#CENT" } != 0 );
      $avg = sprintf( "%.1f", $avg ) if ( $avg != 0 );

      my $avg_unq = 0;
      $avg_unq = $classes{ $key . "#CENT_OCC_UNQ" } / $classes{ $key . "#CENT_UNQ" } if ( $classes{ $key . "#CENT_UNQ" } != 0 );
      $avg_unq = sprintf( "%.1f", $avg_unq ) if ( $avg_unq != 0 );

      print ALL $key . "\t"
        . $classes{$key} . "\t"
        . ( $classes{ $key . "#KNN" } / ( $num_centers * $center_knn ) ) . "\t"
        . $classes{ $key . "#CENT" } . "\t" . "$avg" . "\t"
        . $classes{ $key . "#CENT_UNQ" } . "\t"
        . "$avg_unq" . "\t"
        . $classes{ $key . "#CENT_OCC_MAX" } . "\n";
    }
    close(ALL);

    system("column -t $out_name.classes > t");
    system("mv t $out_name.classes");
  }
  open( RCMD, "| $r_path/R --vanilla --slave" );

  print RCMD "lens<-read.table(\"$out_name.lens\")\n";
  print RCMD "pdf(file=\"$out_name.lens.pdf\",width=8,height=6)\n";
  print RCMD "bp_lens <- boxplot(t(lens),main=\"seq length distribution per center, centers=$num_centers, knn=$center_knn\\n$center_file_name\",xlab=\"CENTER\",ylab=\"seq length\")\n";

  print RCMD "lens_all<-read.table(\"$out_name.lens_all\")\n";
  print RCMD "bp_lens <- boxplot(as.vector(t(lens_all)[,1]),main=\"all seq length distribution, \n$center_file_name\",xlab=\"all seqs\",ylab=\"seq length\")\n";

  # print RCMD "print(as.vector(t(lens_all)[,1]))\n";
  print RCMD "bp_lens <- boxplot(lens_all)\n";
  if ($eval_class) {
    print RCMD "occs<-read.table(\"$out_name.classes\")\n";
    print RCMD "print(as.vector((occs)))\n";
    print RCMD "barplot(names.arg=occs[,1],occs[,3],space=0.2,xlab=\"CLASS\",ylab=\"class density\",main=\"file avg class occ, centers=$num_centers, knn=$center_knn\\n$center_file_name\")\n";
    print RCMD
"barplot(names.arg=occs[,1],occs[,4],space=0.2,xlab=\"CLASS\",ylab=\"#centers\",main=\"center major classes (with first non class 0), centers=$num_centers, knn=$center_knn\\n$center_file_name\")\n";
    print RCMD
"barplot(names.arg=occs[,1],occs[,5],space=0.2,xlab=\"CLASS\",ylab=\"class occ\",main=\"avg class occ for center major class, centers=$num_centers, knn=$center_knn\\n$center_file_name\")\n";
    print RCMD
"barplot(names.arg=occs[,1],occs[,6],space=0.2,xlab=\"CLASS\",ylab=\"#centers\",main=\"center unique major class (without class 0), centers=$num_centers, knn=$center_knn\\n$center_file_name\")\n";
    print RCMD
"barplot(names.arg=occs[,1],occs[,7],space=0.2,xlab=\"CLASS\",ylab=\"class occ\",main=\"avg class occ for center unique major class, centers=$num_centers, knn=$center_knn\\n$center_file_name\")\n";
    print RCMD "barplot(names.arg=occs[,1],occs[,8],space=0.2,xlab=\"CLASS\",ylab=\"class occ\",main=\"max class occ in all centers\\n$center_file_name\")\n";
  }
  print RCMD "dev.off()\n";
  close(RCMD);
}

sub evalSVECTOR {
  my $dir = $_[0];
  my $idx = $_[1];
  my $out = $_[2];

  opendir( DIR, "$dir" );
  my @files = readdir(DIR);
  closedir(DIR);
  my @centers = ();

  foreach my $file (@files) {

    if ( $file =~ /(\d+)\.centers_qual/ ) {

      #print $file. "\n";
      my @lines = readpipe("cat $dir/$file");
      chomp(@lines);
      @lines = map { $_ . " round $1" } @lines;
      push( @centers, @lines );
    }
  }

  # open( IN, "$dir/$idx.centers_qual" );
  #my @centers = <IN>;
  #close(IN);
  #chomp(@centers);

  my @pvals = ();
  my @occ   = ();

  foreach my $idx ( 0 .. $#centers ) {
    my @ent = split( " ", $centers[$idx] );
    $centers[$idx] = \@ent;
    $centers[$idx]->[13] = 0 if ( $centers[$idx]->[13] eq "NaN" );
  }

  @centers = sort { $a->[13] <=> $b->[13] } @centers;

  my @cutoffs = ( 1, 1e-01, 1e-02, 1e-03, 1e-04 );

  open( OUT, ">$out" );

  foreach my $round ( 0 .. $idx ) {

    my @centers_round;
    if ( $round != 0 ) {
      @centers_round = grep { $_->[15] == $round } @centers;
    } else {
      @centers_round = @centers;
    }

    foreach my $co (@cutoffs) {
      my $count = 0;

      #print "new cutoff: $co\n";
      my $occ = 0;
      foreach my $cent (@centers_round) {
        if ( $cent->[13] <= $co ) {
          $count++;
          $occ += $cent->[5];
        }

        #print $cent->[13]." $co\n";

      }

      my $occ_final = $count ? $occ / $count : 0;
      $occ_final = sprintf( "%.1f", $occ_final );
      print "svector_qual round $round pval_cutoff $co centers_rel " . ( $count / scalar(@centers) ) . " centers_abs $count best_class_occurences $occ_final\n";
      print OUT "svector_qual round $round pval_cutoff $co centers_rel " . ( $count / scalar(@centers) ) . " centers_abs $count best_class_occurences $occ_final\n";
    }
  }
  close(OUT);
  makeCol("$out");
}

sub writeAln {
  my $aln_ref = $_[0];
  my $file    = $_[1];
  my $struct  = $_[2];
  my $header  = $_[3];

  my $max_len = 0;
  foreach my $id ( grep { $_ !~ /CONS_STRUCT/ } keys %{$aln_ref} ) {
    $max_len = length($id) if ( length($id) > $max_len );
  }

  open( OUT, ">$file" ) or die "Cannot open aln file $file for writing! Exit...\n\n";

  print OUT "CLUSTAL W --- ";
  if ($header) {
    print OUT "$header\n\n\n";
  } elsif ( exists $aln_ref->{CONS_STRUCT_INFO} ) {
    print OUT $aln_ref->{CONS_STRUCT_INFO} . "\n\n\n";
  } else {
    print OUT "\n\n\n";
  }

  foreach my $seq ( grep { $_ !~ /CONS_STRUCT/ } sort keys %{$aln_ref} ) {
    print OUT sprintf( "%-" . $max_len . "s", $seq ) . "   " . join( "", @{ $aln_ref->{$seq} } ) . "\n";
  }
  if ( $struct && exists $aln_ref->{CONS_STRUCT} ) {
    print OUT sprintf( "%-" . $max_len . "s", " " ) . "   " . join( "", @{ $aln_ref->{CONS_STRUCT} } ) . "\n";
  }
  close(OUT);
}

sub readAln {
  my $aln_file = $_[0];

  open( IN, "$aln_file" ) or die "Cannot open aln file $aln_file for reading! Exit...\n\n";

  my %aln = ();

  while ( my $line = <IN> ) {

    chomp($line);

    next if ( !$line || $line =~ /^\s+$/ );
    if ( $line =~ /^CLUSTAL W -+\s*(.*)\s*$/ ) {
      my $info = $1;
      $info = " " if ( !$info );

      $aln{CONS_STRUCT_INFO} = $info;
      next;
    }

    my @ent = split( " ", $line );

    if ( @ent == 2 ) {
      $ent[1] =~ s/[-\._~]/-/g;    ## count -._~ as similar gap "-" symbol
      $aln{ $ent[0] } .= $ent[1];  #[ split( "", $ent[1] ) ];
    } elsif ( @ent == 1 ) {
      if ( $ent[0] !~ /[\w\*\#]/ && $ent[0] =~ /[\(\)\<\>]/ ) {
        $aln{CONS_STRUCT} .= $ent[0];    #[ split( "", $ent[0] ) ];
      }

    } else {
      print $line. "\n";
      die "Wrong alignment format in $aln_file! Exit...\n\n";
    }
  }

  foreach my $key ( grep { $_ !~ /CONS_STRUCT_INFO/ } keys %aln ) {
    $aln{$key} = [ split( "", $aln{$key} ) ];
  }

  $aln{CONS_STRUCT_BP} = getBasePairs( $aln{CONS_STRUCT} ) if ( exists $aln{CONS_STRUCT} );

  return \%aln;
}

sub getBasePairs {
  my $struct_aref = $_[0];

  my @bp_ret = ();
  map { push( @bp_ret, -1 ) } @{$struct_aref};
  my @bp_stack = ();

  foreach my $col ( 0 .. $#{$struct_aref} ) {

    if ( $struct_aref->[$col] eq "(" ) {
      push( @bp_stack, $col );
    }

    if ( $struct_aref->[$col] eq ")" ) {
      die "Error! Structure is not consistent! Exit...\n\n" if ( !@bp_stack );
      my $open = pop(@bp_stack);
      $bp_ret[$open] = $col;
      $bp_ret[$col]  = $open;
    }

    if ( $struct_aref->[$col] =~ /[^\(\)\-\.\~]/ ) {
      print "WARNING: unknown symbol '" . $struct_aref->[$col] . "' in structure string!\n";
    }
  }

  die "Error! Structure is not consistent! Exit...\n\n" if (@bp_stack);

  return \@bp_ret;
}

sub mlocarna_center {
  my $fasta    = $_[0];
  my $dir      = $_[1];
  my $dpDir    = $_[2];
  my $use_locP = $_[3];

  my $OPTS_locarna_p_model = $CONFIG{OPTS_locarna_p_model};
  my $OPTS_locarna_model   = $CONFIG{OPTS_locarna_model};
  my $loc_path             = $CONFIG{PATH_LOCARNA};

  #  my $vrna_path            = $CONFIG{PATH_VRNA};

  my $loc_pp_dir = "$dir/input";
  system("mkdir -p $loc_pp_dir");

  my @fa = GraphClust::read_fasta_file($fasta);
  foreach my $key ( keys %{ $fa[0] } ) {
    system("ln -f -s $dpDir/$key $loc_pp_dir/$key") if ( -e "$dpDir/$key" && !$use_locP );
  }

  if ($use_locP) {

    ## check mlocarna internal plfold call opts, like prob-cutoff etc.
    ## really use plfold? prob-cutoff is not set in mlocarna and default is 0.01
    system_call( "$loc_path/mlocarna $OPTS_locarna_p_model --verbose --probabilistic --tgtdir $dir $fasta > $dir/locarnaP.out 2>&1 ", 1 );
    my $rel_cmd = "$loc_path/reliability-profile.pl -v -fit-once-on --structure-weight=1  --fit-penalty 0.01 --beta 200 --show-sw --out=$dir/results/result.aln.rel_plot.pdf $dir";
    print $rel_cmd. "\n";
    my @rel = readpipe($rel_cmd);
    open( REL, ">$dir/results/result.aln.rel_signal" );
    print REL @rel[ 0 .. 2 ];
    close(REL);

  } else {

   #system("mlocarna --sequ-local true --tgtdir $dir $fasta $in_mlocarna_opts");
    system_call( "$loc_path/mlocarna $OPTS_locarna_model --verbose --skip-pp --tgtdir $dir $fasta > $dir/locarna.out 2>&1 ", 1 );
  }

  ## alifold result to get consensus structure string for infernal and some nice pictures
#  my $tmp_dir = "$dir/alifold_$$";
#  my $currDir = getcwd;
#  mkdir($tmp_dir);
#  chdir($tmp_dir);
#  my @call_alifold = readpipe("$vrna_path/RNAalifold -r --noLP --color --aln $dir/results/result.aln 2>/dev/null");
#  my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
#  chomp($aln_cons_struct);
#  $aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
#  $aln_cons_struct = $1;
#  open( CONS, ">$dir/results/result.alifold" );
#  print CONS $call_alifold[0];
#  print CONS "$aln_cons_struct\n";
#  print CONS $2;
#  close(CONS);
#  system("mv alirna.ps $dir/results/result.alirna.ps");
#  system("mv aln.ps $dir/results/result.aln.ps");
#  chdir($currDir);
#  system("rm -R $tmp_dir");

  ## cleanup
  system("rm -r -f $dir/input $dir/intermediates $dir/probs $dir/results/single_reliabilities $dir/results/ali* $dir/alifold.out");

}

sub evalCLUSTER {
  my $dir = $_[0];    ## dir with subtree qual/score files
  my $idx = $_[1]; ## idx of round to generate file wth subtree scores in addition
  my $out      = $_[2];    ## outprefix of output files
  my $evaluate = $_[3];    ## evaluate mode

  opendir( DIR, "$dir" );
  my @files = readdir(DIR);
  closedir(DIR);

  my @clusters = ();
  my @scores   = ();
  foreach my $file (@files) {

    if ( $evaluate && $file =~ /bestClusters\.qual\.(\d+)\.\d+/ ) {

      #print $file. "\n";
      my $line = readpipe("head -n 1 $dir/$file");
      chomp $line;
      if ($line) {
        $line .= " round $1";
        push( @clusters, $line );
      } else {
        $line .= "center_id: 1 class_idx: 0 O: 0 C: 0 S: -1 D: 0 [1] 1 round $1";
        push( @clusters, $line );
      }
    }

    if ( $file =~ /bestClusters\.scores\.\d+\.\d+/ ) {

      my @lines = readpipe("cat $dir/$file");
      chomp(@lines);
      push( @scores, @lines );
    }

  }

  ## evaluate cluster quality, pvalue
  if ($evaluate) {
    foreach my $cidx ( 0 .. $#clusters ) {
      my @ent = split( " ", $clusters[$cidx] );
      $clusters[$cidx] = \@ent;
    }

    @clusters = sort { $a->[13] <=> $b->[13] } @clusters;

    my @cutoffs = ( 1, 0.2, 0.1, 0.01 );

    open( OUT, ">$out" );

    foreach my $round ( 0 .. $idx ) {
      my @clusters_round;
      if ( $round != 0 ) {
        @clusters_round = grep { $_->[15] == $round } @clusters;
      } else {
        @clusters_round = @clusters;
      }

      foreach my $co (@cutoffs) {
        my $count = 0;

        #print "new cutoff: $co\n";
        my $perf     = 0;
        my $perf_all = 0;
        foreach my $cent (@clusters_round) {
          if ( $cent->[13] <= $co ) {
            $count++;
            $perf++ if ( $cent->[5] == $cent->[9] );
          }
          $perf_all++ if ( $cent->[5] == $cent->[9] );

          # print $cent->[13]." $co\n";

        }
        my $centers_rel = 0;
        if (@clusters_round) {
          $centers_rel = $count / scalar(@clusters_round);
          $centers_rel = sprintf( "%.3f", $centers_rel );
        }

#print "center_qual round $round pval_cutoff $co centers_rel " . $centers_rel . " centers_abs $count perfect $perf perfect_abs $perf_all\n";
        print OUT "center_qual round $round pval_cutoff $co centers_rel " . $centers_rel . " centers_abs $count perfect $perf perfect_abs $perf_all\n";

      }    ## foreach $co cutoff
      print OUT "\n";
    }    ## foreach round
    close(OUT);

    #system("column -t $out > $out.t");
    #system("mv $out.t $out");
  }

  ## eval subtree scores, do also in non evaluate mode
  my @perfectNodes = ();
  my @noiseNodes   = ();
  my @bestNodes    = ();

  foreach my $line (@scores) {

    # print $line. "\n";
    $line =~ /rank\s+(\S+)\s+/;
    my $rank = $1;
    $line =~ /QUAL_ABS\s+(\S+)/;
    my $qual = $1;

    ## perfect???
    my $perfect = 0;
    $perfect = 1 if ( $evaluate && $qual eq "1.00" );

    #print $line. "\n";

    push( @perfectNodes, [ split( " ", $line ) ] ) if ( $evaluate && $perfect == 1 );
    push( @noiseNodes, [ split( " ", $line ) ] ) if ( $perfect == 0 );
    push( @bestNodes,  [ split( " ", $line ) ] ) if ( $rank == 1 );
  }

  # my $sort_col = 9; ## svm score -> obsolete
  my $sort_col = 45;    ## col after 'SC_SIMPLE'

  @noiseNodes   = sort { $b->[$sort_col] <=> $a->[$sort_col] } @noiseNodes;
  @perfectNodes = sort { $a->[$sort_col] <=> $b->[$sort_col] } @perfectNodes;
  @bestNodes    = sort { $b->[$sort_col] <=> $a->[$sort_col] } @bestNodes;

  open( ALL, ">$out.scores.all" );
  print ALL "TOP RANK NODES:\n";
  print ALL map { join( " ", @{$_} ) . "\n" } @bestNodes;

  if ($evaluate) {
    print ALL "TRUE NODES:\n";
    print ALL map { join( " ", @{$_} ) . "\n" } @perfectNodes;
  }

  if ($evaluate) {
    print ALL "FALSE NODES:\n";
  } else {
    print ALL "ALL NODES:\n";
  }

  print ALL map { join( " ", @{$_} ) . "\n" } @noiseNodes;
  close(ALL);

  #system("column -t $out.scores.all > $out.t");
  #system("mv $out.t $out.scores.all");

  if ($evaluate) {
    open( SVM,    ">$out.scores.SVM" );
    open( SVMABS, ">$out.scores.SVM_abs" );

    ## just to get col names
    my @t;
    map { push( @t, 2 ) } 1 .. 44; ## 44 is to simulate a vector of the expected length for getFeatures
    my @ent = getFeatures( \@t );
    print SVM "QUAL  ";
    print SVM ( join( "   ", sort grep { $_ !~ /QUAL/ } keys %{ $ent[2] } ) ) . "\n";

    print SVMABS "QUAL  ";
    print SVMABS ( join( "   ", sort grep { $_ !~ /QUAL/ } keys %{ $ent[2] } ) ) . "\n";

    ## get now real feature values
    map { @ent = getFeatures($_); print SVM $ent[0] . "\n" } @perfectNodes;
    map { @ent = getFeatures($_); print SVM $ent[0] . "\n" } @noiseNodes;

    map { @ent = getFeatures($_); print SVMABS $ent[1] . "\n" } @perfectNodes;
    map { @ent = getFeatures($_); print SVMABS $ent[1] . "\n" } @noiseNodes;

    close(SVM);
    close(SVMABS);
  }

  open( SUB, ">$out.scores.round$idx" );
  print SUB "TOP RANK NODES:\n";
  print SUB map { join( " ", @{$_} ) . "\n" } grep { $_->[1] =~ /$idx\.\d+/ } @bestNodes;

  if ($evaluate) {
    print SUB "TRUE NODES:\n";
    print SUB map { join( " ", @{$_} ) . "\n" } grep { $_->[1] =~ /$idx\.\d+/ } @perfectNodes;
  }

  if ($evaluate) {
    print SUB "FALSE NODES:\n";
  } else {
    print SUB "ALL NODES:\n";
  }
  print SUB map { join( " ", @{$_} ) . "\n" } grep { $_->[1] =~ /$idx\.\d+/ } @noiseNodes;
  close(SUB);

  #system("column -t $out.scores.round$idx > $out.t");
  #system("mv $out.t $out.scores.round$idx");

}

## create special blacklist for SVM evaluation
## from $round_Scores_file are $bl_groupsize class chooses which are the
## purest/best ones during that round and will be blacklistes in addition for next round
sub evalSVM_blacklist {
  my $data_dir           = $_[0];
  my $round_scores_file  = $_[1];
  my $black_out_filename = $_[2];

  my $bl_ratio     = 1;
  my $bl_groupsize = 1;

  open( CH, "$data_dir/class.hash" );
  my @cl = <CH>;
  close(CH);
  chomp(@cl);

  my %cl2seq = ();
  my %seq2cl = ();

  foreach my $line (@cl) {

    my @ent = split( " ", $line );

    $seq2cl{ $ent[0] } = $ent[1];

    if ( !exists $cl2seq{ $ent[1] } ) {
      my @t = ( $ent[0] );
      $cl2seq{ $ent[1] } = \@t;
    } else {
      push( @{ $cl2seq{ $ent[1] } }, $ent[0] );
    }
  }

  open( SC, "$round_scores_file" );
  my %subtrees = ();
  while ( my $line = <SC> ) {
    chomp($line);
    next if ( $line !~ /class/ );

    $line =~ /cent\s+(\S+)\s+rank\s+(\S+)\s+.*class\s+(\S+)\s+/;
    my $cent  = $1;
    my $rank  = $2;
    my $class = $3;

    $subtrees{"$cent.$rank"} = $class;
  }

  my %subtree_quals = ();

  foreach my $tree ( sort keys %subtrees ) {
    my $str = $subtrees{$tree};
    my @ent = split( ":", $str );

    my %cl_counts = ();
    map { $cl_counts{$_}++ } @ent;
    map { $cl_counts{$_} = $cl_counts{$_} / @ent } keys %cl_counts;
    $cl_counts{"0"} = 1 / @ent if ( exists $cl_counts{"0"} );

    foreach my $cl ( keys %cl_counts ) {
      if ( !exists $subtree_quals{$cl} ) {
        $subtree_quals{$cl} = [ $cl_counts{$cl} ];
      } else {
        push( @{ $subtree_quals{$cl} }, $cl_counts{$cl} );
      }
    }
  }

  my %class_quals = ();

  foreach my $cl ( keys %subtree_quals ) {
    my $qual_avg = 0;
    map { $qual_avg += $_ } @{ $subtree_quals{$cl} };
    $qual_avg = $qual_avg / @{ $subtree_quals{$cl} };
    $class_quals{$cl} = $qual_avg;
  }

  my @cl_qual_sorted = sort { $class_quals{$b} <=> $class_quals{$a} } grep { $_ != 0 } keys %class_quals;
  push( @cl_qual_sorted, 0 ) if ( exists $class_quals{"0"} );
  open( LIST, ">$black_out_filename.list" );
  map { print LIST $_ . " " . $class_quals{$_} . "\n" } @cl_qual_sorted;
  close(LIST);

  open( BL, ">$black_out_filename" );
  foreach my $cl_count ( 1 .. $bl_groupsize ) {

    ## do not blacklist class 0
    next if ( $cl_qual_sorted[ $cl_count - 1 ] == 0 );

    my @ids = @{ $cl2seq{ $cl_qual_sorted[ $cl_count - 1 ] } };
    @ids = shuffle(@ids);
    @ids = @ids[ 0 .. ( int( @ids * $bl_ratio ) - 1 ) ];
    map { print BL $_ . " " . $seq2cl{$_} . "\n"; } sort @ids;
  }
  close(BL);

}

sub subtreeQual {
  my $class_aref = $_[0];
  my $minSeqs    = $_[1];
  my $maxSeqs    = $_[2];

  my %cl_counts = ();
  my %cl_norm   = ();
  my $alpha     = 0.7;

  map { $cl_counts{$_}++ } @{$class_aref};
  my $numcl = keys %cl_counts;
  $cl_counts{"0"} = 1 if ( exists $cl_counts{"0"} );

  #my $norm_size = getNormSize( scalar(@{$class_aref}), $minSeqs, $maxSeqs );
  map { $cl_norm{$_}   = $cl_counts{$_} / $maxSeqs } keys %cl_counts;
  map { $cl_counts{$_} = $cl_counts{$_} / @{$class_aref} } keys %cl_counts;

  my @cl_quals_abs = map { $cl_counts{$_} } sort { $cl_counts{$b} <=> $cl_counts{$a} } keys %cl_counts;
  my @cl_quals_norm =
    map { $alpha * $cl_counts{$_} + ( 1 - $alpha ) * $cl_norm{$_} }
    sort { ( $alpha * $cl_counts{$b} + ( 1 - $alpha ) * $cl_norm{$b} ) <=> ( $alpha * $cl_counts{$a} + ( 1 - $alpha ) * $cl_norm{$a} ) } keys %cl_counts;

  return ( \@cl_quals_norm, \@cl_quals_abs );
}

## calcs a normalized quality of  #$numclass elements (numclass between min/max)
sub getNormSize {
  my $numClass = $_[0];
  my $min      = $_[1];
  my $max      = $_[2];

  my $ymin = 0.2;
  my $ymax = 1;

  #print "min $min max $max ymin $ymin ymax $ymax \n";
  my $a = ( -$ymin + $ymax ) / ( $max - $min );
  my $b = $ymax - $max * ( ( -$ymin + $ymax ) / ( $max - $min ) );

  my $norm = $a * $numClass + $b;

  return $norm;
}

sub normMatrix {
  my $matrix_file = $_[0];
  my $out_prefix  = $_[1];

  return ( $out_prefix . ".normalized", $out_prefix . ".MODEL.norm" );
}

sub getFeatures {
  my $c = $_[0];

  my $rec = {
    QUAL     => $c->[41],
    QUAL_ABS => $c->[43],
    CENT     => $c->[7],
    SCI      => $c->[15],
    MPI      => $c->[17],
    MFE      => $c->[21],
    DENS     => $c->[23],
    BSUM     => $c->[25],
    N        => $c->[29],
    MEAN     => $c->[35],

    #             SD   => $c->[37],
    ALN => $c->[39],
  };

  $rec->{SCI_LOG} = log( $rec->{SCI} ) if ( defined( $rec->{SCI} ) && $rec->{SCI} > 0 );
  $rec->{SCI_LOG} = log(0.01) if ( defined( $rec->{SCI} ) && $rec->{SCI} == 0 );

  $rec->{MPI_LOG} = log( $rec->{MPI} ) if ( defined( $rec->{MPI} ) );
  $rec->{MFE_LOG} = log( -$rec->{MFE} ) if ( $rec->{MFE} && $rec->{MFE} < 0 );
  $rec->{MFE_LOG} = log(0.01) if ( defined( $rec->{MFE} ) && $rec->{MFE} >= 0 );
  $rec->{N_LOG}    = log( $rec->{N} )    if ( $rec->{N} );
  $rec->{DENS_LOG} = log( $rec->{DENS} ) if ( $rec->{DENS} );
  $rec->{BSUM_LOG} = log( $rec->{BSUM} ) if ( $rec->{BSUM} );
  $rec->{MEAN_LOG} = log( $rec->{MEAN} ) if ( $rec->{MEAN} );
  $rec->{ALN_LOG}  = log( $rec->{ALN} )  if ( $rec->{ALN} );

  map { $rec->{$_} = sprintf( "%.4f", $rec->{$_} ) if ( defined( $rec->{$_} ) ) } keys %{$rec};

  my $line;
  my $line_abs;
  $line     = $rec->{QUAL} . " ";
  $line_abs = $rec->{QUAL_ABS} . " ";
  map { $line .= $rec->{$_} . " " if ( defined( $rec->{$_} ) ) } sort grep { $_ !~ /QUAL/ } keys %{$rec};
  map { $line_abs .= $rec->{$_} . " " if ( defined( $rec->{$_} ) ) } sort grep { $_ !~ /QUAL/ } keys %{$rec};

  $line     .= " # QUAL_ABS " . $rec->{QUAL_ABS};
  $line_abs .= " # QUAL_ABS " . $rec->{QUAL_ABS};

  return ( $line, $line_abs, $rec );
}

## classify all feature vectors in $vectors by model $model_f
sub classifySubtree {
  my $vectors = $_[0];
  my $model_f = $_[1];
  my $dir     = $_[2];

  open( NORM, "$model_f.norm" );
  my @norm = <NORM>;
  chomp(@norm);
  close(NORM);

  my @mean = split( " ", $norm[0] );
  my @sd   = split( " ", $norm[1] );
  print "mean: " . join( ":", @mean ) . "\n";
  open( OUT, ">$dir/data.classify" );
  foreach my $vec ( @{$vectors} ) {

    my @data = @{$vec};

#die "I expect 9 features (+1 info col)! Found ".@data."! Exit...\n\n" if (@data != 10);

    my $out_line = "";
    $out_line .= $data[0] . " ";
    my $i = 1;
    map { my $val = ( $_ - $mean[$i] ) / $sd[$i]; $out_line .= "$i:$val "; $i++ } @data[ 1 .. $#mean ];
    print OUT $out_line . "\n";
  }
  close(OUT);

  system("$BIN_DIR/svm_classify $dir/data.classify $model_f $dir/data.predictions");

  open( IN, "$dir/data.predictions" );
  my @pred = <IN>;
  chomp(@pred);
  close(IN);

  map { $pred[$_] = sprintf( "%.3f", $pred[$_] ) } 0 .. $#pred;

  return ( \@pred );
}

sub makeCol {
  my $col_file = $_[0];

#  my $tmp_dir = "_col_tmp_" . $$;
#  mkdir($tmp_dir);
#  system("column -t $col_file > $tmp_dir/t");
#  system("mv $tmp_dir/t $col_file");
#  system("rm -R -f $tmp_dir");

  ## read content and spit into columns
  open(my $in_file,$col_file);
  my @cont = ();
  my $max_col = 0;
  while (my $line = <$in_file>){
    chomp($line);
    my @t = split(" ",$line);
    $max_col = scalar(@t) if ($max_col < scalar(@t));
    push(@cont,\@t);
  }
  close($in_file);


  ## get max length of all columns
  my @len = ();
  map {push(@len,0)}1..$max_col;

  foreach my $ci (0..$max_col-1){
    foreach my $l (@cont){
      if (@{$l} >= $ci+1){
        my $clen = length($l->[$ci]);
        $len[$ci] = $clen if ($len[$ci]<$clen);
      }
    }
  }

  ## print out cols with fixed width
  my $col_spc = 2;

  open(my $out,">$col_file.col_$$");

  foreach my $l (@cont){
    map{printf $out ("%-*s",$len[$_]+$col_spc,$l->[$_]) } 0..@{$l}-1;
    print $out "\n";
  }
  close($out);

  system("mv $col_file.col_$$ $col_file");
}

sub read_partition {
  my $part_file = $_[0];

  open( IN, "$part_file" ) or die "Cannot open partition file $part_file! Exit...\n\n";

  my @part = ();

  while ( my $line = <IN> ) {

    chomp($line);
    my @ent = split( " ", $line );
    push( @part, \@ent );
  }

  close(IN);

  return \@part;
}

sub getRNAzScores {
  my $aln_file  = $_[0];
  my $rnaz_path = $_[1];
  my $rnaz_opts = $_[2];

  #$rnaz_path = "RNAz" if ( !$rnaz_path || !-e $rnaz_path );
  #$rnaz_path .= "/RNAz" if ( $rnaz_path !~ /RNAz$/ );

  print "RNAz path used: $rnaz_path\n";
  my @rnaz_scores = ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

  $rnaz_opts = "" if ( !$rnaz_opts );

  #die "Cannot find RNAz $rnaz_path! Exit...\n" if ( !-e $rnaz_path );
  my $call = readpipe("$rnaz_path $rnaz_opts $aln_file");

  if ( $call =~ /Mean pairwise identity:\s+(\S+)/ ) {
    $rnaz_scores[0] = $1;
    $call =~ /Shannon entropy:\s+(\S+)/;
    $rnaz_scores[1] = $1;
    $call =~ /G+C content:\s+(\S+)/;
    $rnaz_scores[2] = $1;
    $call =~ /Mean single sequence MFE:\s+(\S+)/;
    $rnaz_scores[3] = $1;
    $call =~ /Consensus MFE:\s+(\S+)/;
    $rnaz_scores[4] = $1;
    $call =~ /Energy contribution:\s+(\S+)/;
    $rnaz_scores[5] = $1;
    $call =~ /Covariance contribution:\s+(\S+)/;
    $rnaz_scores[6] = $1;
    $call =~ /Combinations\/Pair:\s+(\S+)/;
    $rnaz_scores[7] = $1;
    $call =~ /Mean z-score:\s+(\S+)/;
    $rnaz_scores[8] = $1;
    $call =~ /Structure conservation index:\s+(\S+)/;
    $rnaz_scores[9] = $1;
    $call =~ /SVM RNA-class probability:\s+(\S+)/;
    $rnaz_scores[10] = $1;
  } else {
    die "Error! Unexpected RNAz output:\n$call\n";
  }

  return \@rnaz_scores;
}

## read fragments from file
## sort according to SEQID, Sorting is important for fragment-overlap function!!!
sub read_fragments {
  my $frag_file = $_[0];

  my @frags = ();
  open( IN, $frag_file ) or die "Cannot find fragment file $frag_file. Exit...\n\n";
  while ( my $line = <IN> ) {
    chomp $line;

    my @ent = split( " ", $line );
    die "Wrong fragment hash format! Exit...\n\n" if ( @ent < 2 );

    my $frag = str2frag( $ent[1] );
    $frag->{VALUE} = $ent[0];

    push( @frags, $frag );
  }
  close(IN);

  @frags = sort { $a->{SEQID} cmp $b->{SEQID} } @frags;

  return \@frags;
}

sub list2frags {
  my $frag_list = $_[0];

  my @frags = ();

  foreach my $fr ( @{$frag_list} ) {
    push( @frags, str2frag($fr) );
  }

  @frags = sort { $a->{SEQID} cmp $b->{SEQID} } @frags;

  return \@frags;
}

sub str2frag {
  my $str = $_[0];

  ## str: frag-key (eg. "SEQ1#2#13#+")
  $str =~ s/\s+//g;
  my @key = split( "#", $str );
  die "Wrong fragment key format! see: $str - Expected 'SEQ1#5#15#+' Exit...\n\n" if ( @key != 4 );

  my $frag = {};
  $frag->{SEQID}  = $key[0];
  $frag->{START}  = $key[1];
  $frag->{STOP}   = $key[2];
  $frag->{STRAND} = $key[3];
  $frag->{KEY}    = $str;

  return $frag;
}

sub fragment_overlap {
  my $fh1           = $_[0];
  my $fh2           = $_[1];
  my $ol_cutoff     = $_[2];
  my $ignore_strand = $_[3];

  my @overlaps = ();

  @{$fh1} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$fh1};
  @{$fh2} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$fh2};

  my $f2_curr_start_idx = 0;
  my $f2_curr_end_idx   = 0;
  my @f2_keys           = ();

 #    map { print "f1 order: $_ " . $fh1->[$_]->{SEQID} . "\n"; } 0 .. $#{$fh1};
 #    map { print "f2 order: $_ " . $fh2->[$_]->{SEQID} . "\n"; } 0 .. $#{$fh2};

  ## check overlaps in more efficient way than all-vs-all
  ## achive this by special sorting of frags
  foreach my $f1 ( 0 .. $#{$fh1} ) {

#     print "$f1 $f2_curr_start_idx f1:" . $fh1->[$f1]->{KEY} ." ". $fh2->[$f2_curr_start_idx]->{KEY}." $f2_curr_start_idx ".@{$fh2}."\n";

    ## check if current (implies all) fragments in fh2 are "greater" that current f1
    next if ( $f2_curr_start_idx >= @{$fh2} || $fh1->[$f1]->{SEQID} lt $fh2->[$f2_curr_start_idx]->{SEQID} );

    ## get all fragments on the same seq in set $fh2
    if ( $fh1->[$f1]->{SEQID} ge $fh2->[$f2_curr_start_idx]->{SEQID} ) {

      # print "range: $f2_curr_start_idx .. $f2_curr_end_idx\n";
      # my $oldend = $f2_curr_end_idx;

      ## set new start-index in $fh2 for correct SEQID
      ## stop at index with equal SEQID or "lower" SEQID (than $f1 SEQID is not in $fh2 )
      while ( $f2_curr_start_idx < @{$fh2} && $fh1->[$f1]->{SEQID} ne $fh2->[$f2_curr_start_idx]->{SEQID} && $fh1->[$f1]->{SEQID} gt $fh2->[$f2_curr_start_idx]->{SEQID} ) { $f2_curr_start_idx++ }

      $f2_curr_end_idx = $f2_curr_start_idx;
      if ( $f2_curr_start_idx >= @{$fh2} || $fh1->[$f1]->{SEQID} lt $fh2->[$f2_curr_start_idx]->{SEQID} ) {
        next;
      }

      ## set new end-index in $fh2 for current SEQID
      while ( $f2_curr_end_idx < @{$fh2} - 1 && $fh1->[$f1]->{SEQID} eq $fh2->[ $f2_curr_end_idx + 1 ]->{SEQID} ) { $f2_curr_end_idx++ }

      ## sort frags according to startpos
      @f2_keys = sort { $fh2->[$a]->{START} <=> $fh2->[$b]->{START} } $f2_curr_start_idx .. $f2_curr_end_idx;

#map {print "skip $_ ".$fh2->[$_]->{SEQID}."\n";} $oldend+1 .. $f2_curr_start_idx-1;
# print "range: $f2_curr_start_idx .. $f2_curr_end_idx\n";
    }

    ## check all pairs, sorting allows to stop if f2 frag starts after frag f1
    foreach my $f2 (@f2_keys) {

#     print "f1:" . $fh1->[$f1]->{SEQID} . " f2:" . $fh2->[$f2]->{SEQID} . "\n";

      next if ( $fh1->[$f1] == $fh2->[$f2] ); ## check if we compare the same object (pointer)

      ## done: if add_reverse=true then cmsearch_scan_reverse should be false, otherwise double hits
      ## wrong because we scon only on original input fasta,
      ## not on additional reverse seqs which were added as fragments if wanted
      ## check for same strand, be careful with cmsearch hits on reverse strand
      next if ( !$ignore_strand && $fh1->[$f1]->{STRAND} ne $fh2->[$f2]->{STRAND} );

      my $start1 = $fh1->[$f1]->{START};
      my $start2 = $fh2->[$f2]->{START};
      my $end1   = $fh1->[$f1]->{STOP};
      my $end2   = $fh2->[$f2]->{STOP};

      next if ( $end2 < $start1 );    ## f2 frag ist before f1
      last if ( $end1 < $start2 ); ## f2 is after f1, due to sorting no overlap is possible anymore

      my $len1 = $end1 - $start1 + 1;
      my $len2 = $end2 - $start2 + 1;

      # 1-2-2-1 ; 1-2-1-2; 2-1-2-1; 2-1-1-2
      my $left        = $start2 - $start1;
      my $right       = $end1 - $end2;
      my $overlap_len = ( abs( $len1 + $len2 ) - abs($left) - abs($right) ) / 2;

      ## use shortest frag as reference length
      my $ol_ref_len = min( $len1, $len2 );
      my $overlap_ratio = $overlap_len / $ol_ref_len;
      $overlap_ratio = sprintf( "%.2f", $overlap_ratio );

      next if ( $overlap_ratio < $ol_cutoff );

      #      print "$f1 $f2  $fh1->[$f1]->{KEY} $fh2->[$f2]->{KEY} overlap!\n";
      push( @overlaps, [ $f1, $f2, $overlap_ratio ] );
    }
  }
  return \@overlaps;
}

sub read_CM_tabfile_ext {
  my $tab_file     = $_[0];
  my $min_bitscore = $_[1];
  my $max_eval     = $_[2];
  my $significance = $_[3];
  my $tab_name     = $_[4];

  open( TABRES, "$tab_file" ) or die "$tab_file tab_result file could not be loaded!";

  my @cm_hit_scores = ();

  my $use_e_values = 0;
  my $lineLast     = "";
  while ( my $line = <TABRES> ) {

    $lineLast = $line;
    next if ( $line =~ /^\#/ );    ## ignore comments in cmsearch output
    chomp $line;
    my @line = split( " ", $line );

    next if ( @line != 9 );

    $use_e_values = 1 if ( $line[7] ne "-" );

    my $strand = "+";
    my $hit    = {
      SEQID    => $line[1],
      START    => $line[2],
      STOP     => $line[3],
      BITSCORE => $line[6],
      EVALUE   => $line[7],
      STRAND   => $strand,
      NAME     => $tab_name
    };

    ## swap start for reverse strand hits, we treat all hits based on + strand
    if ( $line[2] > $line[3] ) {
      $strand        = "-";
      $hit->{START}  = $line[3];
      $hit->{STOP}   = $line[2];
      $hit->{STRAND} = $strand;
    }

    $hit->{KEY} = $hit->{SEQID} . "#" . $hit->{START} . "#" . $hit->{STOP} . "#" . $hit->{STRAND};

    push( @cm_hit_scores, $hit );

    # print join(":",@{$hit})."\n";
    #    $last_id = $line[1];
  }    ## while line <TABRES>
  close TABRES;

  @cm_hit_scores = () if ( $lineLast !~ /^\#/ ); ## simple check that cmsearch finished correctly

  return \@cm_hit_scores if ( !@cm_hit_scores );

  ## score sort hits
  my $co = 0;

  if ( !$use_e_values && $significance < 1 ) {

    @cm_hit_scores = sort { $b->{BITSCORE} <=> $a->{BITSCORE} } @cm_hit_scores;

    open( OUT, ">$tab_file.scores" );
    foreach my $key (@cm_hit_scores) {
      print OUT $key->{BITSCORE} . "\n";
    }
    close(OUT);

    if ( @cm_hit_scores > 10 && $significance < 1 ) { ## todo: check which number is ok for evd fitting
      my $Rcall = readpipe("export R_ENVIRON=$BIN_DIR/Renviron_config; $CONFIG{PATH_R}/Rscript $BIN_DIR/matrixSignificance.R $tab_file.scores $significance 2>$tab_file.R_LOG");

      #print $Rcall;
      $co = $1 if ( $Rcall =~ /cutoff=\s+(.*)$/ );

      #print "fit evd with significance=$significance\n";
    }

    #system("rm -f $tab_file.scores");
    $co = max( $co, $min_bitscore );

    print "found bitscore cutoff=$co ";

  } elsif ($use_e_values) {
    $co = $max_eval;
    print "found evalue cutoff=$co ";
  } else {
    $co = $min_bitscore;
    print "found bitscore cutoff=$co ";
  }

  ## filter hits for score
  my @cm_hit_scores_fil = ();
  print "size hits " . scalar(@cm_hit_scores);
  foreach my $hit (@cm_hit_scores) {

    if ( $use_e_values && $hit->{EVALUE} <= $co ) {
      push( @cm_hit_scores_fil, $hit );
    } elsif ( !$use_e_values && ( ( $co > 0 && $hit->{BITSCORE} >= $co ) || ( $co == 0 && $hit->{BITSCORE} > 0 ) ) ) {
      push( @cm_hit_scores_fil, $hit );
    }

  }

  @cm_hit_scores_fil = sort { $a->{SEQID} cmp $b->{SEQID} } @cm_hit_scores_fil;
  print " size hits fil " . scalar(@cm_hit_scores_fil) . " $tab_file\n";
  ## return array with record refs of filtered hits
  return \@cm_hit_scores_fil;
}

sub evaluate_cm_hits {
  my $classHash   = $_[0];
  my $hits        = $_[1];
  my $min_overlap = $_[2];

  ## guarantee sorting for fragment overlap
  @{$hits} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$hits};

  foreach my $hit ( @{$hits} ) {
    $hit->{CLASS}          = 0;
    $hit->{CLASS_KEY}      = "";
    $hit->{CLASS_OL}       = 0;
    $hit->{CLASS_MULTI}    = [];
    $hit->{CLASS_KEY_DUPL} = 0;
  }

  ## get only overlapping elements back
  ## classOLs aref: [ [idx1,idx2,ol],[idx1',idx2',ol'], ... ]
  my $classOLs = GraphClust::fragment_overlap( $hits, $classHash, $min_overlap );

  my %hit_class_counts   = ();
  my %classKeys_noStrand = ();
  foreach my $ol ( @{$classOLs} ) {

    ## $ol->[0] is array-idx of @{$hits}, $ol->[1] is array-idx of @{$classHash},
    ## current overlapping class key has smaller overlap
    ## [class-idx,SEQ1#10#100,0.2]
    push( @{ $hits->[ $ol->[0] ]->{CLASS_MULTI} }, [ $classHash->[ $ol->[1] ]->{VALUE}, $classHash->[ $ol->[1] ]->{KEY}, $ol->[2] ] );

    ## collect class keys without strand
    ## assume that we have not different signals on strands
    ## store for each class key best (based on bitscore) hit idx
    my $tmp_key = $classHash->[ $ol->[1] ]->{KEY};
    $tmp_key =~ s/[+-]$//;
    if ( !exists $classKeys_noStrand{$tmp_key} || !$hits->[ $ol->[0] ]->{BITSCORE} || $hits->[ $classKeys_noStrand{$tmp_key} ]->{BITSCORE} < $hits->[ $ol->[0] ]->{BITSCORE} ) {
      $classKeys_noStrand{$tmp_key} = $ol->[0];
    }

  }

  ## get uniqe class count over all fragments, including class 0
  foreach my $hit ( sort { $a->{SEQID} cmp $b->{SEQID} } @{$hits} ) {
    my @hcl = ();
    map { push( @hcl, $_->[0] ) } @{ $hit->{CLASS_MULTI} };

    if (@hcl) {
      ## we have signals with class != 0
      @hcl = unique( (@hcl) );
    } else {
      ## we do not have any overlapping signals -> hit is class 0
      @hcl = ("0");
    }

    foreach my $cl (@hcl) {
      $hit_class_counts{$cl}++;
    }
  }

  ## do not count class 0 if set in CONFIG
  delete $hit_class_counts{"0"} if ( $CONFIG{evaluate_class0_as_fp} == 0 );

  ## get majority class first, if two class have same number of hits,
  ## the higher class idx gets first, prevents class 0 first in such a case
  my @class_sort = sort { $hit_class_counts{$b} <=> $hit_class_counts{$a} || $b <=> $a } keys %hit_class_counts;

  my $clusClass = 0;    ## use class 0 as default cluster class
  $clusClass = $class_sort[0] if (@class_sort); ## we could have no overlapping frags with classHash frags

  ## add to all hits the majority class
  map { $_->{CLASS_CLUSTER} = $clusClass } @{$hits};

  ## assign final class of each hit, if clusClass hit is present then use that
  ## otherwise use class key with best overlap
  ## better?: assign key with next best majority class instead of best ol-key class
  ## CLASS_MULTI never contains class 0 frags if we use Graphclust signals (from file class.hash.scan)
  foreach my $hit ( @{$hits} ) {
    my %hcl = ();
    map { $hcl{ $_->[0] } = $_ } @{ $hit->{CLASS_MULTI} };
    if ( exists $hcl{$clusClass} ) {
      $hit->{CLASS}     = $clusClass;
      $hit->{CLASS_KEY} = $hcl{$clusClass}->[1];
      $hit->{CLASS_OL}  = $hcl{$clusClass}->[2];
    } elsif ( @{ $hit->{CLASS_MULTI} } ) {
      ## use best ol frag
      @{ $hit->{CLASS_MULTI} } = sort { $b->[2] <=> $a->[2] } @{ $hit->{CLASS_MULTI} };
      $hit->{CLASS}     = $hit->{CLASS_MULTI}->[0]->[0];
      $hit->{CLASS_KEY} = $hit->{CLASS_MULTI}->[0]->[1];
      $hit->{CLASS_OL}  = $hit->{CLASS_MULTI}->[0]->[2];
    }
  }

  ## identify duplicates of the same class key, mark $hit->{CLASS_KEY_DUPL} = 1
  ## those are not counted as TP in e.g. in f-measure

  foreach my $hidx ( 0 .. $#{$hits} ) {

    next if ( $hits->[$hidx]->{CLASS} eq "0" );

    my $tmp_key = $hits->[$hidx]->{CLASS_KEY};
    $tmp_key =~ s/[+-]$//;
    my $best_idx = $classKeys_noStrand{$tmp_key};
    if ( $hidx != $best_idx ) {
      $hits->[$hidx]->{CLASS_KEY_DUPL} = 1;
    }
  }

  #  foreach my $hit ( @{$hits} ) {
  #    my $tmp_key = $hit->{CLASS_KEY};
  #    next if ( $hit->{CLASS} eq "0" );
  #
  #    $tmp_key =~ s/[+-]$//;
  #    my $tmp_key2 = $hits->[ $classKeys_noStrand{$tmp_key} ]->{CLASS_KEY};
  #    $tmp_key2 =~ s/[+-]$//;
  #    if ( exists $classKeys_noStrand{$tmp_key} && $tmp_key eq $tmp_key2 ) {
  #      ## "KEY" is real fragment-key of hit, not of class signal
  #      if ( $hits->[ $classKeys_noStrand{$tmp_key} ]->{KEY} ne $hit->{KEY} ) {
  #        $hit->{CLASS_KEY_DUPL} = 1;
  #      }
  #    }
  #  }

  return ($clusClass);
}

sub getEvalClassInfo {
  my $class_info_dir = $_[0];

  my $class_frags = GraphClust::read_fragments( $class_info_dir . "/class.hash.scan" );
  my $class_names = GraphClust::read_partition( $class_info_dir . "/class.size" );
  my $class_frags_data = GraphClust::read_partition( $class_info_dir . "/class.hash" );

  my %classSize      = ();
  my %classNames     = ();
  my %classFragsData = ();

  #foreach my $classFrag ( @{$class_frags} ) {
  #  $classSize{ $classFrag->{VALUE} }++;
  #}

  foreach my $nm ( @{$class_names} ) {
    $classSize{ $nm->[1] }  = $nm->[0];
    $classNames{ $nm->[1] } = $nm->[4];
    $classSize{"0"} = $nm->[0] if ( $nm->[1] eq "0" );
  }

  foreach my $fr ( @{$class_frags_data} ) {
    $classFragsData{ $fr->[0] } = $fr->[1];
  }

  #  $eval_num_class = keys %classSize;
  ## True negatives depends on number of input seqs, not fragments
  ## this is correct as we mesaure F/R on scan over input set
  #  $eval_num_seqs = @{$class_frags};

  return ( $class_frags, \%classSize, \%classNames, \%classFragsData );
}

########################################
## read_clustalw_alignment($fh)
##
## read multiple alignment in CLUSTALW aln format from filehandle $fh
##
## $fh file handle open for reading
##
## returns tuple of
##     * ref of alignment hash,
##         which  associates names to alignment strings
##     * ref of names list
##         giving information about order (first occurence) of names
##         in the input
##     * flag, whether a clustal header was detected
##
########################################
sub read_clustalw_alignment {
  my ($fh) = @_;

  my %aln;

  my @names = ();    ## keep a list of names for order

  my $line;

  my $clustal_header = 0;
  if ( ( $line = <$fh> ) =~ /^CLUSTAL/ ) {
    $clustal_header = 1;
    $line           = <$fh>;
  }

  do {
    if ( $line =~ /^([^\s]+)\s+(.+)/ ) {
      my $name = $1;
      my $seq  = $2;

      if ( !exists $aln{$name} ) {
        push @names, $name;
      }

      $aln{$name} .= $seq;

    }
  } while ( $line = <$fh> );

  return ( \%aln, \@names, $clustal_header );
}

########################################
## read_clustalw_alnloh($filename)
##
## read multiple alignment in CLUSTALW aln format from file $filename
##
## returns ref of alignment list of hash
##
########################################
sub read_clustalw_alnloh {
  my ($aln_filename) = @_;

  my $fh;

  open( $fh, "$aln_filename" )
    || die "Can not read alignment file $aln_filename\n";
  my ( $aln, $names, $header ) =
    read_clustalw_alignment($fh);
  close $fh;

  if ( !$header ) { die "Missing CLUSTAL header in $aln_filename."; }

  ## generate list of hash alignment
  my @alnloh;

  for my $name (@$names) {
    my $seq = $aln->{$name};
    my $entry = { name => $name, seq => $seq };
    push @alnloh, $entry;
  }

  return \@alnloh;
}

############################################################
##
## Parse stockholm format.
##
## @param fh file handle
##
## @returns ref to hash with 'Gx entries' and array of sequences (key
## "seqs") The latter array contains hashs with entries "name", "seq"
##
## @note In case of failure, returns 0 and sets $errormsg to a
## printable error description
##
############################################################
sub read_stockholm {
  my $fh = shift;

  my %rec = ();    # result record

  if ( <$fh> !~ /^# STOCKHOLM 1.0$/ ) {
    die("Input not in Stockholm 1.0 format. Missing header \"# STOCKHOLM 1.0\".\n");

    #return 0;
  }

  my %seqs  = ();    ## hash of sequences
  my @names = ();    ## list of sequence names

  while ( my $line = <$fh> ) {

    if ( $line =~ /^\/\/$/ ) { ## end reading at // marker
      last;
    } elsif ( $line =~ /^\s+$/ ) {    ## skip empty lines
                                      # skip
    } elsif ( $line =~ /^#=(G[CR])\s(\S*)\s+(.*)$/ ) { ## handle per column and per residue information
      my $key = "$1 $2";
      my $val = $3;
      chomp $val;
      if ( !exists $rec{$key} ) {
        $rec{$key} = "";
      }
      $rec{$key} .= $val;
    } elsif ( $line =~ /^#=(\S\S)\s(\S*)\s+(.*)$/ ) { ## handle other general information (keep newlines)
      my $key = "$1 $2";
      my $val = $3;
      if ( !exists $rec{$key} ) {
        $rec{$key} = $val;
      } else {
        $rec{$key} .= "\n" . $val;
      }
    } elsif ( $line =~ /^(\S+)\s+(.*)/ ) {    ## handle sequence information
      my $key = $1;
      my $seq = $2;
      if ( !exists $seqs{$key} ) {
        push @names, $key;
        $seqs{$key} = "";
      }
      chomp $seq;
      $seqs{$key} .= $seq;
    }
  }

  my @seqs_entries = ();
  for my $name (@names) {
    my %entry = ();
    $entry{"name"} = $name;
    $entry{"seq"}  = $seqs{$name};
    push @seqs_entries, \%entry;
  }

  $rec{"seqs"} = \@seqs_entries;

  return \%rec;
}

############################################################
##
## Write clustal format.
##
## @param fh file handle
## @param aln ref to stockholm record
##
## @note Does not write GR entries!
## @note Break lines of sequences (with hardcoded width see below)
##
############################################################
sub write_clustal {
  my $fh  = shift;
  my $aln = shift;

  my $width = 120;

  print $fh "CLUSTAL W\n\n";

  my @keys = keys %$aln;

  my $len      = length( $aln->{"seqs"}[0]{"seq"} );
  my $startpos = 0;

  my $name_len = 0;
  for my $seq ( @{ $aln->{"seqs"} } ) {
    $name_len = length( $seq->{"name"} ) if ( $name_len < length( $seq->{"name"} ) );
  }

  while ( $startpos < $len ) {
    print $fh "\n";

    for my $seq ( @{ $aln->{"seqs"} } ) {
      my $subseq = substr( $seq->{"seq"}, $startpos, $width );
      $subseq =~ tr/\./-/;
      print $fh sprintf( "%-" . ( $name_len + 4 ) . "s", $seq->{"name"} ), " ", $subseq, "\n";
    }

    $startpos += $width;
  }
}

sub aln2alifold {
  my $aln_file  = $_[0];
  my $tmp_path  = $_[1];
  my $vrna_path = $_[2];

  ## alifold result to get consensus structure string for infernal and some nice pictures
  my $tmp_dir = "$tmp_path/alifold_$$";
  my $currDir = getcwd;
  mkdir($tmp_dir);
  chdir($tmp_dir);
  my @call_alifold = readpipe( "$vrna_path/" . "RNAalifold -r --noLP --color --aln $aln_file 2>/dev/null" );
  my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
  chomp($aln_cons_struct);
  $aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
  $aln_cons_struct = $1;
  open( CONS, ">$aln_file.alifold" );
  print CONS $call_alifold[0];
  print CONS "$aln_cons_struct\n";
  print CONS $2;
  close(CONS);
  system("mv alirna.ps $aln_file.alirna.ps");
  system("mv aln.ps $aln_file.ps");
  chdir($currDir);
  system("rm -R $tmp_dir");
}

sub mergeFrags {
  my $fa_ref = $_[0];
  my $del    = $_[1]; ## aref with ids, overlaps on these are delete from result

  my %to_delete = ();
  map { $to_delete{$_} = 1 } @{$del};

  my %frags = ();

  foreach my $id ( keys %{ $fa_ref->[0] } ) {
    ## get orig key from header
    $fa_ref->[2]->{$id} =~ /\s*(\S+#\d+#\d+#\S)\s*/;
    ## fill SEQID, START, STOP, KEY, STRAND
    $frags{$id}               = GraphClust::str2frag($1);
    $frags{$id}->{FRAGID}     = $id;
    $frags{$id}->{SEQ}        = $fa_ref->[0]->{$id};
    $frags{$id}->{MERGED}     = $frags{$id}->{KEY};
    $frags{$id}->{MERGED_STR} = $id . ":" . $frags{$id}->{KEY};
    $fa_ref->[2]->{$id} =~ /(ORIGID.*)$/;
    $frags{$id}->{ORIGHEAD} = $1;
  }
  my @frag_list = sort { $a->{SEQID} cmp $b->{SEQID} } values(%frags);

  if ( @{$del} ) {
    my @frag_list_del = grep { exists $to_delete{ $_->{FRAGID} } } sort { $a->{SEQID} cmp $b->{SEQID} } values(%frags);
    my $ols = GraphClust::fragment_overlap( \@frag_list_del, \@frag_list, 0.1, 1 );

    ## delete all blacklisted frags and their direct overlapping frags
    map { delete $frags{ $_->{FRAGID} } } @frag_list_del;
    map { delete $frags{ $frag_list[ $_->[1] ]->{FRAGID} } } @{$ols};
  }

  ## map all frags to their orig SEQID
  my %olmap = ();
  @frag_list = sort { $a->{SEQID} cmp $b->{SEQID} } values(%frags);
  map { $olmap{ $_->{SEQID} }->{ $_->{FRAGID} } = $_ } @frag_list;

  my @newFrags = ();

  foreach my $seq ( keys %olmap ) {

    while ( my @fr = sort { $a->{START} <=> $b->{START} } values %{ $olmap{$seq} } ) {

      my $ols = GraphClust::fragment_overlap( [ $fr[0] ], \@fr, 0.1, 0 );
      #print "check: ".$fr[0]->{KEY}." next: ".$fr[1]->{KEY}."\n";

      foreach my $ol ( @{$ols} ) {

        ## print "MERGE: ".$fr[$ol->[0]]->{START}.":".$fr[$ol->[0]]->{STOP}. " - ".$fr[$ol->[1]]->{START}.":".$fr[$ol->[1]]->{STOP}. " - ";
        my $old0_stop  = $fr[ $ol->[0] ]->{STOP};
        my $old0_len   = length( $fr[ $ol->[0] ]->{SEQ} );

        ## extend frag $ol->[0] with $ol->[1]
        $fr[ $ol->[0] ]->{START} = min( $fr[ $ol->[0] ]->{START}, $fr[ $ol->[1] ]->{START} );
        $fr[ $ol->[0] ]->{STOP} = max( $fr[ $ol->[0] ]->{STOP}, $fr[ $ol->[1] ]->{STOP} );
        $fr[ $ol->[0] ]->{KEY} = $fr[ $ol->[1] ]->{SEQID} . "#" . $fr[ $ol->[0] ]->{START} . "#" . $fr[ $ol->[0] ]->{STOP} . "#" . $fr[ $ol->[0] ]->{STRAND};

        ## print $fr[$ol->[0]]->{START}.":".$fr[$ol->[0]]->{STOP}." ".$fr[$ol->[0]]->{KEY}."\n";
        delete $olmap{$seq}->{ $fr[ $ol->[1] ]->{FRAGID} };

        ## strand aware merge
        ## frags are sorted on forward strand
        ## and we merge always from left to right on forward strand
        ## i.e. on reverse strand we need to put $ol->[1] before $ol->[0], on + strand we put $ol->[1] after $ol->[0]
        if ( $fr[ $ol->[0] ]->{STRAND} eq "+" ) {
          $fr[ $ol->[0] ]->{SEQ} = substr( $fr[ $ol->[0] ]->{SEQ}, 0, $fr[ $ol->[1] ]->{START} - $fr[ $ol->[0] ]->{START} );
          $fr[ $ol->[0] ]->{SEQ} .= $fr[ $ol->[1] ]->{SEQ};
        } else {
          $fr[ $ol->[0] ]->{SEQ} = substr( $fr[ $ol->[0] ]->{SEQ}, $old0_stop - $fr[ $ol->[1] ]->{START} + 1, $old0_len - ( $old0_stop - $fr[ $ol->[1] ]->{START} ) - 1 );
          $fr[ $ol->[0] ]->{SEQ} = $fr[ $ol->[1] ]->{SEQ} . $fr[ $ol->[0] ]->{SEQ};
        }
        $fr[ $ol->[0] ]->{MERGED_STR} .= ";" . $fr[ $ol->[1] ]->{FRAGID} . ":" . $fr[ $ol->[1] ]->{MERGED};
      } ## foreach ol $fr[0] with other frags on same seq

      ## if there are no overlaps anymore we have final merged frag
      if ( !@{$ols} ) {
        delete $olmap{$seq}->{ $fr[0]->{FRAGID} };
        $fr[0]->{FRAGID} = $fr[0]->{FRAGID} . ".1";
        push( @newFrags, $fr[0] );
      }

    } ## while we have frags on seq

  }

  my %newSeqs  = ();
  my %newHeads = ();

  foreach my $nf (@newFrags) {
    $newSeqs{ $nf->{FRAGID} } = $nf->{SEQ};
    $newHeads{ $nf->{FRAGID} } = $nf->{KEY} . " MERGED " . $nf->{MERGED_STR} . " " . $nf->{ORIGHEAD};
  }

  my @ord = map { $_->{FRAGID} } @newFrags;

  return [ \%newSeqs, \@ord, \%newHeads ];
}

1;

#    foreach my $f2 (sort { ( $a =~ /[^#]+#(\d+)/ )[0] <=> ( $b =~ /[^#]+#(\d+)/ )[0] || $a cmp $b } @f2_keys) {

#sub read_CM_tabfile {
#  my $tab_file = $_[0];
#  my $unique   = $_[1];
#
#  open( TABRES, "<$tab_file" ) or die "$tab_file tab_result file could not be loaded!";
#  my %scores  = ();
#  my $idx     = 0;
#  my $last_id = "";
#  while (<TABRES>) {
#
#    next if (/^\#/);
#
#    chomp;
#    my @line = split;
#    $idx++;
#    ## save only best hit per seq
#    next if ( $unique && $line[1] eq $last_id );
#
#    ## [0:seq-id,1:start,2:stop,3:bit-score,4:e-value]
#    $scores{$idx} = [ $line[1], $line[2], $line[3], $line[6], $line[7] ];
#    $last_id = $line[1];
#  }
#  close TABRES;
#  return \%scores;
#}
