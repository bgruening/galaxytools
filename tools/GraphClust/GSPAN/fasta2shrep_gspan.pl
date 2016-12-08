##!/usr/bin/perl
#use feature ':5.10';
use File::Basename;
use lib dirname($0);    # search skript directory for modules
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use Cwd qw(getcwd abs_path);
use File::Temp qw(tempdir);
use File::Copy;
use POSIX qw(ceil);


=head1 NAME


fasta2shrep_gspan.pl -fasta mysequences.fasta -wins "50,100,150,200" -shift 5 -M 8

=head1 SYNOPSIS

Options:

		HELP
        -help   brief help message
        -man    full documentation

        COMPULSORY
        -fasta	<STRING> e.g. "sequence.fasta"
        		All sequences in fasta format.

       	OPTIONS
        -wins		[INTEGER] e.g. "50,100,200"
        		A list of window sizes to use.
        		If none are given (empty string ''), then the entire sequence is
        		taken with no windows. Each window > 1 required!
        -shift		<INTEGER> e.g. 20
        		The shift of the window, relative to the window size given in
        		percent. So you give which percent of the window size shall be
        		used for the shift. Of course the shift is rounded down to the
        		nearest whole number.
        		Example 20 % of a window 150 would result in a step size of 30 nt.
        		It is a relative parameter, as you can give different window sizes.
        		If you do not give this parameter there is a default shift of 1 nt.
        -cue		Crop unpaired ends.
        		If you give this flag, then the unpaired ends of each
        		single structure are ignored. E.g. the structure
        		...(((...))).. becomes just (((...)))
        -stack		Adds stacking information to graphs. This adds an additional
        		vertex (type P) for each pair of stacked base-pairs and four edges
        		(type p) from each of the involved bases to the new vertex.
        -e		<FLOAT> e.g. 5.0
        		Energy range in kcal/mol (RNAshapes)
        		Use only one of -e and -c!
        -c		<INTEGER> e.g. 10
        		Relative energy range, i.e. percentage (%) of MFE energy (RNAshapes)
        		Use only one of -e and -c!
        -t		<INTEGER> [1-5] e.g. 3 OR "3=0,4=100,5=200"
        		The shape type (RNAshapes). Default is 3.
        		With the list format, the shape level can be changed for different window length
        		"4=100" means that shape level 4 is used from length 100nt (window length)
        		The first given length has to be 0! Not continuous given levels are allowed!
        -M		<INTEGER> e.g. 10
        		Max number of shreps that should be taken per window.
        -u 		Ignore unstable structures (RNAshapes).
        		This option filters out closed structures with positive free energy.
        -r		Calculate structure probabilities for shreps (RNAshapes)
        -i		<INT> e.g. 10
        		Turn on structure sampling and gives number of sampling iterations.
        		Default no sampling (i=0)
        -sample-len	<INT> e.g. 100
        		Only in sampling mode: Sampling is only used for seqs/windows >= given length,
        		Default: sample all lengths (0), if -i > 0
        -q		Turn on shape probabilities for RNAshapes, no sampling mode allowed
        -Tp		<FLOAT> e.g 0.001
        		Filter cutoff for shape probabilities, applied before -M filter!
        -seq-graph-win	add for each window a graph which contains no structure
        -seq-graph-t	add for each 't #' a graph which contains no structure
        -seq-graph-alph change the alphabet of unstructured graphs
        -annotate		<STRING> annotation.tab
        				A file with annotations to be added as abstract graphs
        				on the sequence leven (if given) and on the structure
        				(SHREP) level. The format is has the following TAB-delimited
        				columns: SEQID, START, END, NAMESPACE#LABEL.
        				Labels with the same name-space and SEQID form connected
        				components, which is a sequence of label vertices ordered
        				by the START position in the sequence.
        -abstr			Add abstract structure graphs to the single shrep graph
        				instances.
        -nostr			Calculate no structures, only add sequence information,
        				if this is given, then -seq-graph-win AND/OR -seq-graph-t
        				are required.
        -match-shape    <SHAPE>
                all seqs/windows will be constraint folded into that shape via
                RNAshapes (if structure is given in another way this struct will be kept),
                if this shape is not possible within given energy range, produce a
                specific t graph with only one vertex 'X'. By this the instance
                becomes very unsimilar to all other graphs (for knn)
        -vp     enable graph computation with viewpoints:
                svmsgdnspdk will center on those nucleotides that are given
                via capital letters and ignore those given as lowercase letters
        -tmp		<STRING> e.g. "/scratch/1/sita/tmp"
        		A directory for writing temporary files
        -o		<STRING> e.g. "ProjectX/MySequences/GSPAN/"
        		Output directory for gspan files containing graphs.
        -group		<INTEGER> e.g. 5
                        Combine/group that number of input seqs into 1 gspan file
                        output name is then '<INT>.group.gspan.bz2'

        -stdout         send graphs to stdout instead of writing to files
        -ignore-header  don't write fasta id part after first space to gspan
        -debug          additional debug output


        DEFAULT VALUES
        -wins	""
        -shift	1 nt
        -c		10
        -t		3
        -M		0 # selects all shreps
        -tmp    "/var/tmp/fasta2shrep"
        -o		"CURRENT_DIR/GSPAN_Outputs/"



=head1 DESCRIPTION

=cut

###############################################################################


###############################################################################
# PARSE COMMAND LINE OPTIONS
###############################################################################

# command line options
my ( $i_help, $i_man, $i_debug, $i_fas, $i_wins, $i_shift, $i_crop_unpaired_ends,
  $i_r, $i_e, $i_c, $i_t, $i_u, $i_M, $i_o, $i_sge, $i_jobid, $i_tmp,
  $i_ignore_seq_header, $i_stacks, $i_stdout, $i_q, $i_T, $i_i,
  $i_sample_min_length, $i_sge_logDir, $i_sge_errDir, $i_groupsize, $i_annotate,
  $i_abstr, $i_no_structure, $i_vp, $i_matchShape );

my ( $i_add_seq_graph_win, $i_add_seq_graph_t, $i_change_seq_graph_alph );

my $options = GetOptions(
  "help"            => \$i_help,
  "man"             => \$i_man,
  "debug"           => \$i_debug,
  "fasta=s"         => \$i_fas,
  "wins=s"          => \$i_wins,
  "shift=f"         => \$i_shift,
  "cue"             => \$i_crop_unpaired_ends,
  "stack"           => \$i_stacks,
  "r"               => \$i_r,
  "e=f"             => \$i_e,
  "c=i"             => \$i_c,
  "t=s"             => \$i_t,
  "u"               => \$i_u,
  "M=i"             => \$i_M,
  "tmp=s"           => \$i_tmp,
  "o=s"             => \$i_o,
  "i=i"             => \$i_i,
  "sample-len=i"    => \$i_sample_min_length,
  "q"               => \$i_q,
  "Tp=f"            => \$i_T,
  "seq-graph-win"   => \$i_add_seq_graph_win,
  "seq-graph-t"     => \$i_add_seq_graph_t,
  "seq-graph-alph"  => \$i_change_seq_graph_alph,

  "task-id=i"       => \$i_jobid,
  "ignore-header"   => \$i_ignore_seq_header,
  "stdout"          => \$i_stdout,

  "group=i"         => \$i_groupsize,
  "annotate=s"      => \$i_annotate,
  "abstr"           => \$i_abstr,
  "nostr"           => \$i_no_structure,
  "vp"              => \$i_vp,
  "match-shape=s"   => \$i_matchShape,
);

pod2usage( -exitstatus => 1, -verbose => 1 ) if $i_help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $i_man;
($options) or pod2usage(2);

# check compulsory options
($i_fas) or pod2usage("Error: the option -fasta is compulsory!\n");
( -e $i_fas ) or pod2usage("Error: no such file - $i_fas!\n");
$i_fas = abs_path($i_fas);

# check other options and set default values
pod2usage("Error: either -e OR -c can be given, but NOT both!\n") if ( $i_e && $i_c );
( $i_e or $i_c ) or $i_c = 10;    # set -c 10, if neither -e or -c are given
($i_M) or $i_M = 0;    # max number of shreps is 0 (=means take all computed)
($i_i) or $i_i = 0;    # default no sampling else sampling iterations
($i_sample_min_length) or $i_sample_min_length = 0;
pod2usage("\nError: use --sample-len <INT> only with -i <1..INT> !\n") if ( not $i_i and $i_sample_min_length );
($i_T) or $i_T = 0;
($i_q) or $i_q = 0;
pod2usage("\nError: Sampling (-i) not possible with shape probabilities (-q)!\n") if ( $i_i and $i_q );
($i_add_seq_graph_win)     or $i_add_seq_graph_win     = 0;
($i_add_seq_graph_t)       or $i_add_seq_graph_t       = 0;
($i_change_seq_graph_alph) or $i_change_seq_graph_alph = 0;
($i_no_structure)          or $i_no_structure          = 0;

if ($i_change_seq_graph_alph) {
  ($i_add_seq_graph_t) or ($i_add_seq_graph_win) or pod2usage( "Error: " .
"When giving the parameter -seq-graph-alph, then either -seq-graph-t or -seq-graph-win" .
      " must also be given!\n" );
}
( -e $i_annotate ) or pod2usage("Error: no such file - $i_annotate!\n")
  if ($i_annotate);

$i_add_seq_graph_t = 1 if ($i_no_structure);

($i_t) or $i_t = 3;    # default abstraction type is 3
my $change_shape_level = 0;
my @level_lens = ( -1, -1, -1, -1, -1 ); ## array_idx-1=shape level, value=start length of this level

if ( $i_t !~ /^\d+$/ ) {

  my @t_minlens = split( ",", $i_t );    ## -t "3=100,4=200"

  foreach my $idx ( 1 .. @t_minlens ) {

    my $level = $t_minlens[ $idx - 1 ];    ##  $level = "3=100"

    die "$level Wrong -t format! Example: -t 3=0,4=100,5=200\n" if ( $level !~ /^\d+\=\d+$/ );
    my @lev_len = split( "=", $level );    ##  $level = "3=100"

    die "Wrong -t format! First level given needs to be length 0! Example: -t 3=0,4=100,5=200\n" if ( $idx == 1 && $lev_len[1] != 0 );
    die "Wrong -t format! Only level 1-5 allowed! Example: -t 3=0,4=100,5=200\n" if ( $lev_len[0] < 1 or $lev_len[0] > 5 );
    die "Wrong -t format! Length >= 0 expected! Example: -t 3=0,4=100,5=200\n" if ( $lev_len[1] < 0 );
    $change_shape_level = 1;
    $level_lens[ $lev_len[0] - 1 ] = $lev_len[1];
  }

  ($i_debug) and print STDERR "i_t = $i_t - change $change_shape_level -  shape level lengths#" . join( ":", @level_lens ) . "#" . join( ":", @t_minlens ) . "#\n";
}

## checks for match-shape
if ($i_matchShape){
  die "Please provide correct match shape string like '[]'! Exit...\n\n"
    if ($i_matchShape !~ /^[\[\]_]+$/);
  ## RNAshapes prodices anyway only 1 structure, no suboptimal structs in match-shape folding
  $i_M = 1;
}

my $CURRDIR = getcwd;

# set up tmp directory
# default tmp is /var/tmp, usually not NFS mounted!
( defined $i_tmp ) or $i_tmp = '/var/tmp/';
my $tmp_template = 'fasta2shrep-XXXXXX';

# CLEANUP => 1 : automatically delete at exit
$i_tmp = tempdir( $tmp_template, DIR => $i_tmp, CLEANUP => 1 );

# create GSPAN directory when not printing to stdout
if ( not $i_stdout ) {
  if ($i_o) {
    ( -e $i_o ) or system("mkdir -p $i_o");
  } else {
    system("mkdir -p GSPAN_Outputs");
    $i_o = $CURRDIR . "/GSPAN_Outputs/";
  }
}

###############################################################################
# GLOBAL VARIABLES
###############################################################################

my $rnashapes_loc = "./RNAshapes";

if ( !$rnashapes_loc || !-e $rnashapes_loc ) {
  my $loc = `which RNAshapes`;
  chomp($loc);
  print " loc = $loc\n";

  die "\n 1. Cannot find RNAshapes binary! Exit...\n\n" if ( !$loc );
  die "\n 2. Cannot find RNAshapes binary! Exit...\n\n" if ( !-e $loc );
  $rnashapes_loc = $loc;

}

my @WINDOWS = ();
@WINDOWS = split( ",", $i_wins ) if ($i_wins);
my $globalFolding;
$globalFolding = 1 unless @WINDOWS;
my $CURRUSER = getlogin;
my $SEQNO    = 0;          # used to id the sequences
my $GSPANNO  = 0;          # used for gspan filenames
# minimum length of sequence that will work with RNAshapes
# sequences with one or two nucleotides will be restricted to sequence-only graphs
# this ensures that no sequences are skipped and external info kept synchronized
my $GSPAN_SEQ_MINLEN = 3;

# name spaces
my $ABSTRUCT = "AS";
###############################################################################
# EXECUTION CODE
###############################################################################

# read fasta file into hash
my ( $headers_aref, $sequences_aref, $metainfo_aref ) = read_fasta_with_nonunique_headers_meta($i_fas);



my @used_seq_headers;
my @used_seqs;
my @used_meta;
my $group_idx;

if ($i_jobid) {
  ## just process the one sequence as given by the jobid number
  my $used_grouping = 1;   ## if no group is given, make 1 seq per job = group=1
  $used_grouping = $i_groupsize if ($i_groupsize);
  my $st  = ( $i_jobid - 1 ) * $used_grouping + 1;
  my $end = $st + $used_grouping - 1;
  $end = @{$sequences_aref} if ( $end > @{$sequences_aref} );

  foreach my $idx ( $st .. $end ) {
    push( @used_seq_headers, $headers_aref->[ $idx - 1 ] );
    push( @used_seqs,        $sequences_aref->[ $idx - 1 ] );
    push( @used_meta,        $metainfo_aref->[ $idx - 1 ] );
  }

  $group_idx = $i_jobid;

  #($i_debug) and print STDERR "st $st end $end gr $group_idx job $i_jobid gs $i_groupsize\n";

  $GSPANNO = $i_jobid-1 if (!$i_groupsize);

} else {
  ## process all sequences at once
  @used_seq_headers = @{$headers_aref};
  @used_seqs        = @{$sequences_aref};
  @used_meta        = @{$metainfo_aref};
}

my $out;
my $gspanfile;
my $out_no_match_shape;

if ($i_matchShape && !$i_stdout && !$i_groupsize){
  open($out_no_match_shape,">$i_o/fasta2shrep.no_match");
}

# for each sequence in the fasta file
while ( my $seq = shift @used_seqs ) {
  my $tmp_header = shift @used_seq_headers;
  my $tmp_meta   = shift @used_meta;

  my ( $seq_id, $seq_header ) = ( $tmp_header =~ /(\S+)\s*([\S*\s*]*)/ );
  $i_ignore_seq_header and $seq_header = '';

  my $seq_fasta = generate_single_fasta_from_sequence_X( $seq_id, $seq );
  my $seq_len = length($seq);

  # only print sequence graphs for sequences below this threshold
  my $no_structure_override = ($seq_len < $GSPAN_SEQ_MINLEN) ? 1 : 0;

  $GSPANNO++;

  # set outstream for gspan output to correct file/STDOUT
  if ($i_stdout) {
    $out = \*STDOUT;
  } elsif ( !$i_groupsize ) {
    $gspanfile = $i_tmp . '/' . $GSPANNO . '.gspan';
    open( $out, "| bzip2 -f > $gspanfile.bz2" );
  } elsif ( ( $GSPANNO - 1 ) % $i_groupsize == 0 ) {

    if ( $GSPANNO > 1 ) {
      close($out);
      system("mv $gspanfile.bz2 $i_o/$group_idx.group.gspan.bz2");
      if ($out_no_match_shape){
        close($out_no_match_shape);
        system("mv $gspanfile.no_match $i_o/$group_idx.group.gspan.no_match");
      }
    }

    if ( !$i_jobid ) {
      $group_idx = int( ( $GSPANNO - 1 ) / $i_groupsize ) + 1;
    }

    $gspanfile = "$i_tmp/$group_idx.group.gspan";
    open( $out, "| bzip2 -f > $gspanfile.bz2" );
    open( $out_no_match_shape, ">$gspanfile.no_match" ) if ($out_no_match_shape);
  }

  ## do not use folding windows in special cases
  if ( $globalFolding || $i_no_structure || $no_structure_override) {
    @WINDOWS = ();
    push( @WINDOWS, $seq_len );
  }

  ##check win sizes for global folding
  my @WINDOWS_used = sort { $a <=> $b } @WINDOWS;
  foreach my $w_idx ( 0 .. ( @WINDOWS_used - 1 ) ) {
    if ( $WINDOWS_used[$w_idx] >= $seq_len && $WINDOWS_used[$w_idx] > 1 ) {
      if ( $w_idx < $#WINDOWS_used ) {
        @WINDOWS_used = @WINDOWS_used[ 0 .. $w_idx ];
        last;
      }
    }
  }

  ## use seq graph only if no shape folding wanted
  $i_add_seq_graph_t = 1 if ($i_no_structure || $no_structure_override);

  ## no shape info in graphheader if we have fixed structure (according LocaRNA handling)
  ## use tags #FS or #S for provided structure
  my $graph_header;
  my @struct_meta = grep { $_ =~ /#FS/ || $_ =~ /#S/ } keys %{$tmp_meta};

  if (@struct_meta) {
    $graph_header = getGraphHeader( $seq_id, $seq_header, \@WINDOWS, $i_shift, $i_e, $i_c, $i_t, $i_u, $i_r, $i_M, $i_crop_unpaired_ends, $i_i, $i_sample_min_length, $i_q, $i_T, $seq_len, 1 );
  } else {
    $graph_header = getGraphHeader( $seq_id, $seq_header, \@WINDOWS, $i_shift, $i_e, $i_c, $i_t, $i_u, $i_r, $i_M, $i_crop_unpaired_ends, $i_i, $i_sample_min_length, $i_q, $i_T, $seq_len, $i_no_structure || $no_structure_override );
  }

  print $out $graph_header;

  my $gi = 0;    # current graph index

  ## add graph with no structure at all depending on $i_add_seq_graph_t
  if ($i_add_seq_graph_t) {

    ($gi) = convertSeqWindow( $seq, $seq_len, 1, $gi, $graph_header, $out, $i_annotate, $seq );

  }

  ## encode fixed structure only if provided and structures wanted in general
  if ( @struct_meta && !$i_no_structure && !$no_structure_override) {

    my $struct_meta = "#FS";
    $struct_meta = "#S" if ( !exists $tmp_meta->{$struct_meta} );

    my $seq_shrep = [ $tmp_meta->{$struct_meta}, "ENERGY", "0.00", "SHAPE", $struct_meta ];
    $gi = convertShapeWindow( [$seq_shrep], $seq, $seq_len, 1, $gi, $out,
      $graph_header, $i_annotate, $i_abstr, $i_crop_unpaired_ends, $i_stacks, $seq );
    @WINDOWS_used = ();    ## no shape folding if we have a fixed structure

  }

  ## ignore RNAshapes folding if wanted, but do correct file move afterwards
  ## (just "next" does not work due to output)
  @WINDOWS_used = () if ($i_no_structure or $no_structure_override);

  #for each window size in list
  foreach my $win_size (@WINDOWS_used) {

    # calculate shift size from percentage
    my $curr_shift = 1;
    if ($i_shift) {
      $curr_shift = ( $i_shift / 100 ) * $win_size;
      $curr_shift = int($curr_shift);                 #round down
      $curr_shift = 1 unless ($curr_shift);           # just in case it is 0
    }
    ($i_debug) and print STDERR "winsize: $win_size curr_shift: $curr_shift\n";
    ($i_debug) and print STDERR "\nNext: $seq_id\t winsize:$win_size \n";

    # choose current shape level, depending on $i_t
    my $curr_t = 0;
    if ($change_shape_level) {
      for ( my $i = 0 ; $i < @level_lens ; $i++ ) {
        $curr_t = $i + 1 if ( $level_lens[$i] != -1 && ( $level_lens[$i] <= $win_size ) );
      }
      ($i_debug) and print STDERR "$win_size curr type $curr_t\n";
    } else {
      $curr_t = $i_t;
    }

    my $rnashapesoutput_fh;

    # call RNAshapes and write to $rnashapesoutput_fh
    $rnashapesoutput_fh = call_RNAshapes( $seq_fasta, $rnashapes_loc, $win_size,
      $curr_shift, $i_e, $i_c, $curr_t, $i_u, $i_r, $i_q, $i_T, $i_i, $i_sample_min_length, $seq_len, $i_matchShape );

    # read RNAshapes output from $rnashapesoutput_fh and write subgraph
    # to gspan file
    my $gi_old = $gi;
    $gi = convert_RNAshapes_output( $rnashapesoutput_fh, $gi, $i_M, $out, $graph_header,
      $win_size, $seq_len, $curr_t, $i_annotate, $i_abstr, $i_crop_unpaired_ends, $i_stacks, $seq );

    ## no (match) shape found at all for this seq
    if ($gi == $gi_old+1){
      $gi = convertSeqWindow( "X", 1, 1, $gi, $graph_header, $out, $i_annotate, $seq );
      print $out_no_match_shape $seq_id."\n" if (!$i_stdout && $out_no_match_shape);
    }
  }    ## foreach WINDOW_used

  if ( !$i_stdout && !$i_groupsize ) {
    close($out);
    move "$gspanfile.bz2", "$i_o/$GSPANNO.gspan.bz2";
  }

  system("rm $seq_fasta");

}    ## while @used_seqs

if ($i_groupsize) {
  close($out);
  system("mv $gspanfile.bz2 $i_o/$group_idx.group.gspan.bz2");
  close($out_no_match_shape) if ($out_no_match_shape);
  move "$gspanfile.no_match", "$i_o/$group_idx.group.gspan.no_match" if ($out_no_match_shape);
} elsif ( $out_no_match_shape && !$i_stdout && !$i_groupsize ) {
  close($out_no_match_shape);
}

###############################################################################
# METHODS
###############################################################################

############################################################################
# Generates fasta file for a single sequence (one-lined-fasta). This
# fasta is stored in the temp directory and should be deleted at the end.
# Input:
# seq_id : the sequence ID
# seq : the sequence
#
# Output:
# The fasta file name
############################################################################
sub generate_single_fasta_from_sequence_X {
  my ( $seq_id, $seq ) = @_;

  $seq = uc($seq);
  $seq =~ tr/T/U/;
  $seq =~ s/[^AUCGN]/N/g;

  my $outfas = $i_tmp . "/seq_" . $SEQNO++ . ".fasta";
  my $host   = readpipe("hostname");
  open( FAS, ">$outfas" ) or die "$host Cannot open file $outfas! Exit...\n\n";
  print FAS ">$seq_id\n$seq";
  close(FAS);

  return $outfas;
}

############################################################################
# RNAshapes is called with the given input or default parameters.
# Input:
# seq_fasta : the sequence fasta file
# rnashapes_location : the location of the installation files for RNAshapes
# win_size : the current window size
# shift : the input parameter -shift
# e : the input parameter -e
# c : the input parameter -c
# t : the input parameter -t
# u : the input parameter -u
# r : the input parameter -r
#
# Output: none
############################################################################
sub call_RNAshapes {
  my ( $seq_fasta, $rnashapes_location, $win_size, $shift, $e, $c, $t, $u, $r, $q, $T, $i, $sample_length, $seqLen, $matchShape ) = @_;
  my $FUNCTION = "call_RNAshapes in fasta2shrep_gspan.pl";

  ($seq_fasta) or die("INPUT ERROR in $FUNCTION: the fasta file is compulsory!\n");
  ($rnashapes_location) or die( "INPUT ERROR in $FUNCTION: the RNAshapes location" . " is compulsory!\n" );
  die "$rnashapes_location does not exist! Exit...\n\n" if ( !-e $rnashapes_location );
  my $call = $rnashapes_location . " -o 1 ";    # the output format is of type 1
  #$call .= "--mode sample ";
  #$call .= "--outputLowProbFilter 0 ";
  $call .= "-q " if ($q);
  $call .= "-T $T " if ( $q and $T );
  $call .= "-w $win_size ";  ##### before it was just --windowSize
  $call .= "-W $shift " if ($i_shift);  ##### before it was just --windowIncrement

  die("ERROR in $FUNCTION: Give only one of the options -c or -e (RNAshapes)!\n")

    if ( $e && $c );
  $call .= "-e $e " if ($e);  ##### before it was just --relativeDeviation
  $call .= "-c $c " if ($c); ##### before it was just --absoluteDeviation
  $call .= "-t $t " if ($t);		##### before it was just --shapeLevel
  $call .= "-u "    if ( $u and not $i );       ## not possible in sampling mode
  $call .= "-r "    if ($r);  ##### before it was just --structureProbs
  $call .= "-m $matchShape " if ($matchShape);

  ## check is a bit long but : we want to sample if the window is larger than $sample_length and full seq is larger than window or sample_len
  ## necessary to do sampling a large window is given, sample_length is shorter than window, but seq is longer than sample_len
  $call .= "-i $i -A " if ( not $q and $i and ( $win_size >= $sample_length and ( $win_size <= $seqLen or $seqLen >= $sample_length ) ) ); ## -A is to omit samples and print only combined shape probs
  $call .= " < $seq_fasta";

  ($i_debug) and print STDERR "$seqLen $sample_length $win_size $call\n";

  open my $rnashapesoutput, "$call |" or die( "ERROR in $FUNCTION: The following call " . "could not be carried out!\n$call\n" );

   print " \n";
   print " call function = 	$call	";
   print "\n ";

  return $rnashapesoutput;
}

############################################################################
# The output of RNAshapes for one sequence and one window size is read
# and converted into graph format.
# Input:
# rnashapeoutput : filehandle for RNAshapes output in format -o 1
# curr_gi : current graph index (for vertices)
# maxShreps : max number of shreps to convert to graphs
# graph_file_hdl : the output handler for the graph file
# graphHead : the header line for the complete graph (for sequence)
# winSize : current window size in input for RNAshapes
# seqLen  : full input seq length
# used_t  :
# annotate:
# abstr   :
# cue     :
# stacks  :
# orig_seq : the nucleotide sequence as read from fasta
#
# Output:
# The current graph index
############################################################################
sub convert_RNAshapes_output {
  my ( $rnashapesoutput, $curr_gi, $maxShreps, $graph_file_hdl, $graphHead, $winSize, $seqLen, $used_t, $annotate, $abstr, $cue, $stacks, $orig_seq ) = @_;

  ## omit first line in output, contains fasta header line
  my $line       = <$rnashapesoutput>;
  my $win_shreps = [];
  my $win_start;
  my $win_end;
  my $win_seq;
  my $win_shrep_count = 0;
  my $win_sample      = 0;
  my $winHead;
  my $win_globalFolding;
  my $win_size_real;

  # reading RNAshapes output
  while ( $line = <$rnashapesoutput> ) {

    if ( $line =~ /^(\d+)\s+(\d+)$/ ) {
      ## line: "<start>    <end>"

      if ( @{$win_shreps} > 0 ) {

        print $graph_file_hdl $winHead;

        ## remove SHAPE="_" shrep if it exists in $win_shreps
        ## do this only when have a sequence graph already
        if ($i_add_seq_graph_t) {
          my @new_win_shreps = ();
          map { push( @new_win_shreps, $_ ) if ( $_->[ $#{$_} ] ne "_" ) } @{$win_shreps};
          $win_shreps = \@new_win_shreps;
        }

        ## add graph with no structure depending on $i_add_seq_graph_win to win
        if ($i_add_seq_graph_win ) {
          ( $curr_gi ) = convertSeqWindow( $win_seq, $win_size_real,
            $win_start, $curr_gi, $winHead, $graph_file_hdl, $annotate, $orig_seq );
        }

        $curr_gi = convertShapeWindow( $win_shreps, $win_seq, $win_size_real,
          $win_start, $curr_gi, $graph_file_hdl, $winHead, $annotate, $abstr,
          $cue, $stacks, $orig_seq );
      }

      ## set new window params
      $win_shreps      = [];
      $win_start       = $1;
      $win_end         = $2;
      $win_shrep_count = 0;
      $win_size_real   = $win_end - $win_start + 1;

      if ( ($win_size_real) >= $seqLen ) {
        $win_globalFolding = 1;
      } else {
        $win_globalFolding = 0;
      }
      my $win_center = $win_start + ( ( $win_size_real + 1 ) / 2 );

      $winHead = getWindowHeader( $graphHead, $winSize, $win_start, $win_end, $win_center, $win_globalFolding, $win_sample, $used_t );

    } elsif ( $line =~ /^(\S+)$/ ) {
      ## line: "CUUAUGAGUAAGGAAAAUAACGAUUCGGGGUGACGCCCGAAUCCUCACUG"
      $win_seq = uc($1); ## to be 101% sure that we have by default uppercase chars

    } elsif ( $line =~ /^Results for (\d+) iterations:$/ ) {
      ## line: "Results for 10 iterations:"
      $win_sample = $1;

    } elsif ( $line =~ /^Shape\s+\S+\s+not found within energy range.*$/ ) {
      ## line: "Shape [] not found within energy range (-24.75 to -27.50). Try -c or -e to increase range."
      $win_shreps      = [];

    } elsif ( $line =~ /^([\(\)\.]+)\s+\((\S+)\)\s+(\S+)$/ ) {
      ## line:"...((((..(((....)))))))...........(((((......)))))  (-10.10)  [[]][]"
      ## take only $maxShreps shreps per window if set
      next if ( $maxShreps && $win_shrep_count >= $maxShreps );
      push( @{$win_shreps}, [ $1, "ENERGY", $2, "SHAPE", $3 ] );
      $win_shrep_count++;

    } elsif ( $line =~ /^([\(\)\.]+)\s+\((\S+)\)\s+\((\S+)\)\s+(\S+)$/ ) {
      ## line:"((((..((((...((.((.((.....)).)).))...))))..))))...  (-10.60)  (0.7795360)  [[[[[]]]]]"
      ## take only $maxShreps shreps per window if set
      next if ( $maxShreps && $win_shrep_count >= $maxShreps );
      push( @{$win_shreps}, [ $1, "ENERGY", $2, "PROB", $3, "SHAPE", $4 ] );
      $win_shrep_count++;

    } elsif ( $line =~ /^([\(\)\.]+)\s+\((\S+)\)\s+(\S+)\s+(\S+)$/ ) {
      ## line:"((((..((((...((.((.((.....)).)).))...))))..))))...  (-10.60) 0.3000000 [[[[[]]]]]"
      ## take only $maxShreps shreps per window if set
      next if ( $maxShreps && $win_shrep_count >= $maxShreps );
      push( @{$win_shreps}, [ $1, "ENERGY", $2, "SHAPEPROB", $3, "SHAPE", $4 ] );
      $win_shrep_count++;

    } elsif ( $line =~ /^([\(\)\.]+)\s+\((\S+)\)\s+\((\S+)\)\s+(\S+)\s+(\S+)$/ ) {
      ## line:"((((..((((...((.((.((.....)).)).))...))))..))))...  (-10.60)  (0.7795360) 0.3000000 [[[[[]]]]]"
      ## take only $maxShreps shreps per window if set
      next if ( $maxShreps && $win_shrep_count >= $maxShreps );
      push( @{$win_shreps}, [ $1, "ENERGY", $2, "PROB", $3, "SHAPEPROB", $4, "SHAPE", $5 ] );
      $win_shrep_count++;

    } else {
      next if ( $line =~ /^$/ );
      die "Unexpected shape output format!\nline=$line\n\nExit...\n\n";
    }
  }

  ## convert last windows
  if ( @{$win_shreps} > 0 ) {

    print $graph_file_hdl $winHead;

    ## remove SHAPE="_" shrep if it exists in $win_shreps
    ## do this only when have a sequence graph already
    if ($i_add_seq_graph_t) {
      my @new_win_shreps = ();
      map { push( @new_win_shreps, $_ ) if ( $_->[ $#{$_} ] ne "_" ) } @{$win_shreps};
      $win_shreps = \@new_win_shreps;
    }

    ## add graph with no structure depending on $i_add_seq_graph_win to win
    if ($i_add_seq_graph_win) {
      ($curr_gi) = convertSeqWindow( $win_seq, $win_size_real, $win_start,
        $curr_gi, $winHead, $graph_file_hdl, $annotate, $orig_seq );
    }

    $curr_gi = convertShapeWindow( $win_shreps, $win_seq, $win_size_real,
      $win_start, $curr_gi, $graph_file_hdl, $winHead, $annotate, $abstr, $cue,
      $stacks, $orig_seq );
  }

  close($rnashapesoutput);

  return $curr_gi + 1;    # return the gi (graph index) for the next subgraph
}

## Sub to create graph for a complete unstructured sequence
# TODO: document function
sub convertSeqWindow {
  my ( $win_seq, $win_size_real, $win_start, $curr_gi, $winHead, $graph_file_hdl, $annotate, $orig_seq ) = @_;

  my $seq_shrep;

  # use a different alphabet (lowercase) if option seq-graph-alph set
  my $seq_graph_sequence = $i_change_seq_graph_alph ? lc($win_seq) : uc($win_seq);
  $seq_shrep = [ "." x $win_size_real, "ENERGY", "0.00", "SHAPE", "_", "STRUCT", "." x $win_size_real, "SEQ", $seq_graph_sequence ];
  my $backboneGraph_ref = getBackboneGraph( $seq_graph_sequence, $curr_gi, $win_start, 0, ( $win_size_real - 1 ), $orig_seq );

  print $graph_file_hdl getSeqHeader( $winHead, $seq_shrep );
  print $graph_file_hdl join( "\n", @{$backboneGraph_ref} ) . "\n";

  $curr_gi += $win_size_real;

  if ($annotate) {
    ## TODO add annotation graphs
  }

  return ($curr_gi);
}

############################################################################
# The results for one window are converted in this method to GSPAN graphs.
# Input:
# win_shreps_aref : the array ref of shreps (dot-bracket structures) for the
# 					current window
# win_seq : the nucleotide sequence for the current window
# win_size_real : the current (true) window size
# win_start : the starting position of win_seq in the original seq (1-n)
# curr_gi : the current graph index
# graph_file_hdl : the graph file handler
# winHead   : the header line for the current sequence window
# TODO annotate
# abstr : input parameter -abstr
# cue   : input parameter -cue
# stacks : input parameter -stacks
# orig_seq : the nucleotide sequence as read from fasta
#
# Output:
# The current graph index
############################################################################
sub convertShapeWindow {
  my ( $win_shreps_aref, $win_seq, $win_size_real, $win_start, $curr_gi,
    $graph_file_hdl, $winHead, $annotate, $abstr, $cue, $stacks, $orig_seq ) = @_;

  ## generate for each shrep a connected component in gspan file
  foreach my $shrep ( @{$win_shreps_aref} ) {

    # get the current gi as it was at the beginning of this shrep
    my $shrep_gi = $curr_gi;

    # cut off unpaired ends, if option is given
    my $shrep_struct     = $shrep->[0];
    my $crop_index_left  = 0;
    my $crop_index_right = length($shrep_struct) - 1;

    # find croping indices with numbering from 0 to n-1
    if ($cue) {
      $crop_index_left = index( $shrep_struct, "(" );    # find 1st occ of "("
      $crop_index_right = rindex( $shrep_struct, ")" );  # find last occ of ")"

      # if the complete window is unpaired, then don't crop
      if ( $crop_index_left == -1 ) {
        $crop_index_left  = 0;
        $crop_index_right = length($shrep_struct) - 1;
      }
    }

    # create structure graph
    my $backboneGraph_ref = getBackboneGraph( $win_seq, $curr_gi, $win_start, $crop_index_left, $crop_index_right, $orig_seq );
    my $structGraph_ref;

    # add both basic edges, stacks and abstract graph
    if ($abstr) {
      ( $structGraph_ref, $curr_gi ) = getStructPlusAbstractStructGraph( $shrep, $curr_gi, $win_size_real, $cue, $stacks );

      # just add basic structure graph edges, and stacks
    } else {
      ( $structGraph_ref, $curr_gi ) = getStructGraph( $shrep, $curr_gi, $win_size_real, $stacks );
    }

    # additional information for shrep header: sequence and dot-bracket
    my $crop_length = $crop_index_right - $crop_index_left + 1;
    my $shrepheader_struct = substr( $shrep_struct, $crop_index_left, $crop_length );
    my $shrepheader_seq = substr( $win_seq, $crop_index_left, $crop_length );
    push( @{$shrep}, 'STRUCT', $shrepheader_struct, 'SEQ', $shrepheader_seq );

    # print structure graph to file
    print $graph_file_hdl getShapeHeader( $winHead, $shrep );
    print $graph_file_hdl join( "\n", @{$backboneGraph_ref} ) . "\n";

    # don't print empty array; sgdnspdk breaks with empty lines
    if ( @{$structGraph_ref} > 0 ) {
      print $graph_file_hdl join( "\n", @{$structGraph_ref} ) . "\n";

      if ($annotate) {

        # TODO add annotations to shrep structure
        # $win_start, $win_end (1-n)
      }

    }
  }

  return $curr_gi;
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









###########################################################################
# Here the vertices for the nucleotides are generated and also the backbone
# edges, which connect the nucleotides in their correct order (i.e order
# of the given sequence)
# Input:
# win_seq : the sequence for the current window
# win_start : the starting position for the current window
# curr_crop_i_left : the cropping index for the left end of sequence (-cue)
# curr_crop_i_right : the cropping index for the right end of sequence (-cue)
# orig_seq : the nucleotide sequence as read from fasta
#
# Output:
# The lines of the graph in an array reference to be printed to a file with
# one element per line.
############################################################################
sub getBackboneGraph {
  my ( $win_seq, $curr_gi, $win_start, $curr_crop_i_left, $curr_crop_i_right, $orig_seq ) = @_;

  # RNAshapes substitutes T -> U, sequence only does not
  # thus, just transform all t/T -> u/U to remain consistent
  $orig_seq =~ tr /tT/uU/;
  $win_seq  =~ tr /tT/uU/;

  my @seq;
  @seq = split( "", $win_seq );

  # if the vp option is set we need to obtain the original capitalization from
  # $orig_seq and extract the windowed sequence manually
  my @seq_vp;
  if ( defined $i_vp ) {
    my $win_len = length($win_seq);
    my $capitalized_win_seq = substr( $orig_seq, $win_start - 1, $win_len );
    ( uc($capitalized_win_seq) eq uc($win_seq) ) or
      die( "error: windowed sequence generated due to vp option not equal " .
        "to sequence reported by RNASHAPES.\n" .
        "'${capitalized_win_seq}' != '${win_seq}'" );
    @seq_vp = split( "", $capitalized_win_seq );
  }

  my @vert = ();
  my @edg  = ();

  my $curr_abs_pos = $win_start + $curr_crop_i_left;
  $curr_gi += $curr_crop_i_left;

  # set vertice labeled with 'v' or 'V' according to sequence and vp option
  # when vp set,
  # uppercase nucleotides are annotated with 'v'
  # lowercase nucleotides are annotated with 'V'
  my $vertice_label = ( defined $i_vp and ( $seq_vp[$curr_crop_i_left] =~ /[a-z]/ ) ) ? 'V' : 'v';

  # create backbone vertice of first nucleotide
  push( @vert, join( ' ', $vertice_label, $curr_gi, $seq[$curr_crop_i_left], $curr_abs_pos ) );

  foreach my $idx ( ( $curr_crop_i_left + 1 ) .. $curr_crop_i_right ) {
    $curr_abs_pos++;
    $curr_gi++;

    # set vertice label as described above
    $vertice_label = ( defined $i_vp and ( $seq_vp[$idx] =~ /[a-z]/ ) ) ? 'V' : 'v';
    push( @vert, join( ' ', $vertice_label, $curr_gi, $seq[$idx], $curr_abs_pos ) );
    push( @edg, join( ' ', "e", $curr_gi - 1, $curr_gi, '> 1' ) );
  }

  my @ret = ( @vert, @edg );

  return \@ret;
}

###########################################################################
# This method does the same as getStructGraph, but is extended to identify
# the abstract parts of the structure and also add this to the graph
# via abstract relations, see information below.
# We already have the backbone graph, now the base-pair edges need to be
# added to this graph and the abstract graph is added after this. Finally
# the abstract relations, pointing from the abstract level to the basic
# structure level is added. All the while, we keep track of the current
# index graph. Furthermore, in the basic structure graph, we have vertices
# stacks that connects the four stacked base-pairs involved. (symbols P
# for vertex and p for edges). The abstract vertices are labelled according
# to the given name-space, then #, then the abstract structure type, e.g.
# HL for hairpin-loop. The edges between the structure types are denoted
# with n, and the relation vertex is ^R, with ^r going from the abstract
# structure to the relation and from the relation to the basic structure
# is @r.
#
# Input:
# curr_shrep : the current shrep as dot-bracket structure
# curr_gi : the current graph index
# win_size : the window size
# cue : the input parameter -cue (whether to crop unpaired ends)
# stacks : the input parameter -stacks (whether to add stack information)
#
# Output:
# The structure graph information line-by-line as an array reference and
# the current graph index.
#
############################################################################
sub getStructPlusAbstractStructGraph {
  my ( $curr_shrep, $curr_gi, $win_size, $cue, $stacks ) = @_;

  #  $win_size =  length($curr_shrep->[0]);

  my @struct = split( "", $curr_shrep->[0] );

#  print "=================testing: new shrep : gi=$curr_gi ===============================\n" if ($i_debug);

  # OBJECTS and VARIABLES
  # all indices are saved according to the current graph index, $curr_gi,
  # which is at the beginning of the current window, so that the nucleotide
  # index can be inferred using $idx + $curr_gi
  # opening brackets
  my @open_blks = (); # open blocks of consecutive opening brackets (array of arrays)
  my @p_open = ();  # currently open block of the last identified complete block
  my @c_open = ();  # current "incomplete" block of consecutive opening brackets
                    # base-pair objects
  my %bps    = ();  # hash of all base-pairs with stem for a value, if it is at
                    # the BP is opening or closing a stem object
  my @c_bp   = ();  # current BP, just being closed at current index pos (i,j)
  my @p_bp   = ();  # previously closed BP at position before current index
                    # stem objects
  my %stems  = ();  # all stems, key="i:j,k:l", value=stem array, (i,j) outer BP
  my @c_stem = ();  # current incomplete stem object
  my @stmbrks = (); # all stem-breaks in the form of (i,k,l), where i is the remaining
                    # open base and (k,l) is the subsequent closing base-pair

#  my @p_stmbrk	= (); # previous stem-break in the form of (i,k,l), where i is the remaining
#  					  # open base and (k,l) is the subsequent closing base-pair
# loop objects
  my @up    = (); # all unpaired regions (i,j) that are of unkown type
  my @c_up  = (); # current incomplete unpaired object
  my @hls   = (); # all hairpin loop objects
  my @bls   = (); # all bulge-loop objects
  my @ils   = (); # all internal loops
  my @mls   = (); # all multi-loops
  my @c_mls = (); # current incomplete multi-loops! (there can be more than one)
  my @els   = (); # all external loops
  my @c_el  = (); # current incomplete external loop

  # iterate through shrep structure an identify abstract parts
  foreach my $idx ( 0 .. @struct - 1 ) {

    #===================================================
    # current char is the opening bracket of a base-pair
    if ( $struct[$idx] eq "(" ) {

      #===================================================

      # update current window index to current graph index
      $idx += $curr_gi;

      # case: '((' currently there is an open block, extend it
      if (@c_open) {
        push( @c_open, $idx );

# case: '.(' or ')('  or BOS'(' there is no open block, create a new one ".(" or ")("
      } else {

        # begin a new current open block
        push( @c_open, $idx );

        # case: ')('
        if (@c_stem) {

          # CLOSE STEM
          print STDERR "TEST - closed stem:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@c_stem ) ) if ($i_debug);
          close_stem( \@c_stem, \%stems, \@p_bp, \@c_bp, \@p_open, \@stmbrks, \@open_blks, \%bps );

          # case (x(x)(
          if (@open_blks) {

            # EXTEND_ML, if one exists
            if (@c_mls) {
              extend_or_open_nested_ML( $curr_shrep, \@c_mls, \@p_bp, \@up, \@stmbrks, \@open_blks, $idx, $curr_gi );

              # OPEN_ML
            } else {

              # case 'x((x)(' the opening of the ML was a stem-break
              if ( @stmbrks
                && $stmbrks[-1]->[1] == $p_bp[0]
                && $stmbrks[-1]->[2] == $p_bp[1] ) {
                my @c_stmbrk = @{ pop(@stmbrks) };
                my @a        = ( $c_stmbrk[0], -1 );
                my @ml       = ();
                push( @ml, \@a );    # 1st base-pair is in-waiting
                my @a2 = ( $c_stmbrk[1], $c_stmbrk[2] );
                push( @ml, \@a2 );    #2nd BP
                my @a3 = ($idx);
                push( @ml,    \@a3 );    #first of current BP
                push( @c_mls, \@ml );
                print STDERR "TEST - opened ML with SB:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@ml, 1 ) ) if ($i_debug);

                # case 'x(.(x)('
              } elsif (@up) {
                my @tmp_up = @{ pop(@up) };

                # compare adjacent BP to UP to see if they fit
                if ( $tmp_up[1] + 1 == $p_bp[0]
                  && $open_blks[-1]->[-1] == $tmp_up[0] - 1 ) {
                  my @a = ( $open_blks[-1]->[-1], -1 );
                  my @ml = ();
                  push( @ml, \@a );    # 1st BP awaits closing
                  my @newbp = @p_bp;
                  push( @ml, \@newbp );    # 2nd BP
                  my @a2 = ($idx);
                  push( @ml,    \@a2 );    # current open BP, awaits closing
                  push( @c_mls, \@ml );
                  print STDERR "TEST - opened ML with prev UP1:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@ml, 1 ) )

                    if ($i_debug);
                } else {
                  die "ERROR: the base-pairs to not match the " . "previous unpaired region?!\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
                }
              } else {
                die "ERROR: in case ML, but there is no initial" . "part of the ML to be identified\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
              }
            }

            # case: 'x(x)('
          } else {

            # case: '.(x)(' => CLOSE_EL
            if (@c_el) {
              if ( $c_el[-1] == $p_bp[0] ) {
                push( @c_el, $p_bp[1] );
                my @newel = @c_el;
                push( @els, \@newel );
                @c_el = ();
                print STDERR "TEST - closed EL1:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newel ) ) if ($i_debug);
              } else {
                die "ERROR: the previous BP must match the current EL\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
              }
            }

            # in both cases OPEN_EL
            push( @c_el, @p_bp );
            push( @c_el, $idx );
            print STDERR "TEST - opened EL:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@c_el ) ) if ($i_debug);
          }

          # case '.('
        } elsif (@c_up) {

          # case (x(x).( OR 'x(.('
          if (@open_blks) {

            # case 'x(x).(' ML
            if ( @p_bp && ( $p_bp[1] == $c_up[0] - 1 ) ) {

              # case '((x)x(x).(' EXTEND_ML, if one exists
              if (@c_mls) {

                extend_or_open_nested_ML( $curr_shrep, \@c_mls, \@p_bp, \@up, \@stmbrks, \@open_blks, $idx, $curr_gi );

                # case 'x((x).(' or 'x(.(x).(' OPEN_ML
              } else {

                # case 'x((x).(' the opening of the ML was a stem-break
                if ( @stmbrks
                  && $stmbrks[-1]->[1] == $p_bp[0]
                  && $stmbrks[-1]->[2] == $p_bp[1] ) {
                  my @c_stmbrk = @{ pop(@stmbrks) };
                  my @a        = ( $c_stmbrk[0], -1 );
                  my @ml       = ();
                  push( @ml, \@a );    # 1st base-pair is in-waiting
                  my @a2 = ( $c_stmbrk[1], $c_stmbrk[2] );
                  push( @ml, \@a2 );    #2nd BP
                  my @a3 = ($idx);
                  push( @ml,    \@a3 );    #first of current BP
                  push( @c_mls, \@ml );
                  print STDERR "TEST - open ML with SB:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@ml, 1 ) ) if ($i_debug);

                  # case 'x(.(x).(' opening the ML with initial UP region
                } elsif (@up) {

                  # get previous unpaired region
                  my $tmp  = pop(@up);
                  my @p_up = @{$tmp};

          # case: '(.(x).(' compare adjacent BP to UP to see if they create a ML
                  if ( $p_up[1] + 1 == $p_bp[0]
                    && $open_blks[-1]->[-1] == $p_up[0] - 1 ) {
                    my @a = ( $open_blks[-1]->[-1], -1 );
                    my @ml = ();
                    push( @ml, \@a );    # 1st BP awaits closing
                    my @newbp = @p_bp;
                    push( @ml, \@newbp );    # 2nd BP
                    my @a2 = ($idx);
                    push( @ml,    \@a2 );    # current open BP, awaits closing
                    push( @c_mls, \@ml );
                    print STDERR "TEST - opened ML with prev UP2:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@ml, 1 ) )

                      if ($i_debug);
                  }

                  # has to be an ML, but can't find first part
                } else {
                  die "ERROR: in case ML, but there is no initial " . "part of the ML to be identified\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
                }
              }

              # case '(.(
            } else {
              my @newup = @c_up;
              push( @up, \@newup );
              print STDERR "TEST - closed UP:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newup ) ) if ($i_debug);
            }

            # case: 'x(x).(' OR 'EOS.('
          } else {

            # case: '.(x)(' => CLOSE_EL
            if (@c_el) {
              if ( $c_el[-1] == $c_bp[0] ) {
                push( @c_el, $p_bp[1] );
                my @newel = @c_el;
                push( @els, \@newel );
                @c_el = ();
                print STDERR "TEST - closed EL2:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newel ) ) if ($i_debug);
              } else {
                die "ERROR: the previous BP must match the current EL\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
              }
            }

            # in both cases OPEN_EL
            if (@p_bp) {
              push( @c_el, @p_bp );
              push( @c_el, $idx );
            } else {

              # ignore this if -cue is given
              # TODO check the gi index
              if ($cue) {
                @c_el = ();
              } else {
                push( @c_el, ( $curr_gi, $curr_gi ) );
                push( @c_el, $idx );
                print STDERR "TEST - opened EL:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@c_el ) ) if ($i_debug);
              }
            }
          }
          @c_up = ();

        } else {

          # no more cases except opening bracket at beginning of sequence
          # in this case do nothing, as the index has already been added
          die "ERROR: '((', ')(', '.(' have all been covered\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) )
            unless ( $idx == $curr_gi );

        }
      }

      #===================================================
      # current char is the unpaired base
    } elsif ( $struct[$idx] eq "." ) {

      #===================================================

      # update current window index to current graph index
      $idx += $curr_gi;

      # case ').' => CLOSE_STEM
      if (@c_stem) {
        print STDERR "TEST - closed stem:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@c_stem ) ) if ($i_debug);
        close_stem( \@c_stem, \%stems, \@p_bp, \@c_bp, \@p_open, \@stmbrks, \@open_blks, \%bps );
      }

      # case '..' extend UP, @c_up
      if (@c_up) {
        $c_up[1] = $idx;

        # case 'x.' open new UP
      } else {

        # OPEN_UP, @c_up
        @c_up = ( $idx, $idx );

        # case '(.' just come to the end of an open block, push it onto stack
        if (@c_open) {
          my @newopen = @c_open;
          push( @open_blks, \@newopen );
          @c_open = ();

          # case ').' or '.'
          # (have already closed stem and created stem-break and emptied p_open)
        } else {

          # case 'x(x).' or '.' => CLOSE_EL
          if ( !@open_blks ) {

            # case 'x().().'
            if (@c_el) {

              # double check
              if ( $c_el[-1] == $p_bp[0] ) {
                push( @c_el, $p_bp[1] );
                my @newel = @c_el;
                push( @els, \@newel );
                @c_el = ();
                print STDERR "TEST - closed EL3:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newel ) ) if ($i_debug);
              } else {
                die "This case should not occur!\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
              }
            }
          }
        }
      }

      #===================================================
      # current char is the closing bracket of a base-pair
    } elsif ( $struct[$idx] eq ")" ) {

      #===================================================

      # update current window index to current graph index
      $idx += $curr_gi;

      # save previous base-pair
      if (@c_bp) {
        @p_bp = @c_bp;
      }

      # case: '((x))', extend stem
      if (@p_open) {

        # get current base-pair
        @c_bp = ( pop(@p_open), $idx );

        # add it to the base-pair hash
        $bps{"$c_bp[0]:$c_bp[1]"} = "";

        # EXTEND_STEM
        if (@c_stem) {
          $c_stem[0] = $c_bp[0];
          $c_stem[1] = $c_bp[1];
        } else {
          die "ERROR: in case '))' there has to be an open stem object!\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
        }

        # case: '(x(x))' or '(x.)', '()' not allowed
      } else {

        # case: '(x(x))'  CLOSE_STEM close previous stem
        if (@c_stem) {
          print STDERR "TEST - closed stem:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@c_stem ) ) if ($i_debug);
          close_stem( \@c_stem, \%stems, \@p_bp, \@c_bp, \@p_open, \@stmbrks, \@open_blks, \%bps );

          # case '(x(x))' OPEN_STEM
          if (@open_blks) {
            my $tmp = pop(@open_blks);
            @p_open = @{$tmp};

            # get current base-pair
            @c_bp = ( pop(@p_open), $idx );

            # add it to the base-pair hash
            # as the stem is not finished yet, cannot add stem information
            $bps{"$c_bp[0]:$c_bp[1]"} = "";

            # open stem
            $c_stem[0] = $c_bp[0];
            $c_stem[1] = $c_bp[1];
            $c_stem[2] = $c_bp[0];
            $c_stem[3] = $c_bp[1];
          } else {
            die "ERROR: there have to be open blocks to match " . "current closing bracket!\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
          }

          # case: '(.(x))' => CLOSE_BL
          if ( @up && $up[-1]->[0] - 1 == $c_bp[0] && $up[-1]->[1] + 1 == $p_bp[0] ) {
            my @newBL = ( @c_bp, @p_bp );
            push( @bls, \@newBL );
            pop(@up);
            print STDERR "TEST - closed BL:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newBL ) ) if ($i_debug);

            # case: '(x(x).(x))' => CLOSE_ML
          } elsif ( @c_mls && $c_mls[-1]->[-1]->[-1] == $p_bp[0] ) {
            my @newml = @{ pop(@c_mls) };
            push( @{ $newml[-1] }, $p_bp[1] );    # close last BP
            $newml[0]->[1] = $c_bp[1];            # close 1st BP
            push( @mls, \@newml );
            print STDERR "TEST - closed ML1:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newml, 1 ) ) if ($i_debug);
          } else {
            die "ERROR: what is this case?\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
          }

          # case: (x.)
        } else {

          # OPEN_STEM open new stem
          if (@open_blks) {
            @p_open = @{ pop(@open_blks) };

            # get current base-pair
            @c_bp = ( pop(@p_open), $idx );

            # add it to the base-pair hash
            # as the stem is not finished yet, cannot add stem information
            $bps{"$c_bp[0]:$c_bp[1]"} = "";

            # open stem
            $c_stem[0] = $c_bp[0];
            $c_stem[1] = $c_bp[1];
            $c_stem[2] = $c_bp[0];
            $c_stem[3] = $c_bp[1];
          } else {
            die "ERROR: there have to be open blocks to match " . "current closing bracket!$curr_shrep->[0]\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
          }

          # CLOSE_C_UP
          if (@c_up) {

            # case: (.) => close_HL
            if ( $c_up[0] - 1 == $c_bp[0] ) {
              my @newhl = @c_bp;
              push( @hls, \@newhl );
              @c_up = ();
              print STDERR "TEST - closed HL:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newhl ) ) if ($i_debug);

              # case: (x(x).)
            } elsif ( @p_bp && $c_up[0] - 1 == $p_bp[1] ) {

              # case: ((x).) => CLOSE_BL
              if ( @stmbrks
                && $stmbrks[-1]->[0] == $c_bp[0]
                && $stmbrks[-1]->[1] == $p_bp[0]
                && $stmbrks[-1]->[2] == $p_bp[1] ) {
                my @newbl = ( @c_bp, @p_bp );
                push( @bls, \@newbl );
                pop(@stmbrks);
                @c_up = ();
                print STDERR "TEST - closed BL:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newbl ) ) if ($i_debug);

                # case: (x(x).), either CLOSE_ML or CLOSE_IL
              } else {

                # case '(.(x).)' => CLOSE_IL
                if ( @up
                  && ( $up[-1]->[0] - 1 == $c_bp[0] )
                  && ( $up[-1]->[1] + 1 == $p_bp[0] ) ) {
                  my @p_up = @{ pop(@up) };
                  my @newil = ( @c_bp, @p_bp );
                  push( @ils, \@newil );
                  @c_up = ();
                  print STDERR "TEST - closed IL:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newil ) ) if ($i_debug);

                  # case (x(x)x(x).) => CLOSE_ML remember, array of tuples
                } elsif ( @c_mls && $c_mls[-1]->[-1]->[0] == $p_bp[0] ) {
                  my @newml = @{ pop(@c_mls) };
                  push( @{ $newml[-1] }, $p_bp[1] );    # close final BP
                  $newml[0]->[1] = $c_bp[1];            # close 1st BP
                  push( @mls, \@newml );
                  @c_up = ();
                  print STDERR "TEST - closed ML2:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@newml, 1 ) ) if ($i_debug);
                } else {
                  die "There should be no such case!\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
                }
              }

              # there should be no more unknown unpaired regions left!
            } else {
              die "There shouldn't be any unknown unpaired " . "regions left\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
            }

          } else {
            die "ERROR: '()' is not allowed in $curr_shrep->[0]\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) );
          }
        }
      }
    }
  }

  # EOS

  # CLOSE_STEM
  print STDERR "TEST - closed stem EOS:\n" . ( test( $curr_shrep->[0], ( $win_size + $curr_gi ), $curr_gi, \@c_stem ) )
    if ( @c_stem && $i_debug );
  close_stem( \@c_stem, \%stems, \@p_bp, \@c_bp, \@p_open, \@stmbrks, \@open_blks, \%bps ) if (@c_stem);
  die "ERROR: there should not be an open ML!\n" . $curr_shrep . "\n"
    if (@c_mls);

  # case ().() CLOSE_EL
  if (@c_el) {
    if ( $c_el[-1] == $p_bp[0] ) {
      push( @c_el, $p_bp[1] );
      my @newel = @c_el;
      push( @els, \@newel );
      @c_el = ();

      print STDERR "TEST - closed EL4 EOS:\n" . ( test( $curr_shrep->[0], ( $win_size + $curr_gi ), $curr_gi, \@newel ) )
        if ($i_debug);
    } else {
      die "ERROR: the base-pairs in the EL don't match" . $curr_shrep . "\n";
    }
  }

  # case (). CLOSE_EL
  if (@c_up) {
    unless ($cue) {
      if (@p_bp) {
        push( @c_el, @p_bp );
      } else {
        push( @c_el, ( $curr_gi, $curr_gi ) );
      }

      push( @c_el, ( ( $curr_gi + $win_size - 1 ), ( $curr_gi + $win_size - 1 ) ) );
      my @newel = @c_el;
      push( @els, \@newel );
      @c_el = ();
      print STDERR "TEST - closed EL5 EOS:\n" . ( test( $curr_shrep->[0], ( $win_size + $curr_gi ), $curr_gi, \@newel ) )
        if ($i_debug);
      ## check case where cue is given, but window is completely unstructured
    } else {
      unless (@p_bp) {
        my @newel = ( $curr_gi, $curr_gi, ( $curr_gi + $win_size - 1 ), ( $curr_gi + $win_size - 1 ) );
        push( @els, \@newel );
      }
    }
  }

  #  # update current graph index (since all shrep indices have been read)
  $curr_gi += $win_size;

  # ===================================================
  # build graph lines
  my @graph_lines = ();

  # ===================================================

  # graph line objects
  my @edg        = ();
  my @stackgraph = ();
  my @abstrgraph = ();

 #-------------------> stems,stacks
 # build basic graph edges, stacks, and vertices for stems in the abstract graph
  my @tmp = keys %stems;
  foreach my $sk (@tmp) {
    my @stem = @{ $stems{$sk} }; # stem (i,j,k,l)=> (i,j)= (0,1) and (k,l)=(2,3)
    die "ERROR: There are not enough elements in this stem: ", ( join( ",", @stem ) ), "\n" unless ( scalar(@stem) == 4 );

    # infer base-pairs from stem
    @p_bp = ();
    @c_bp = ();
    my @stem_relation_idxs  = ();
    my @stack_relation_idxs = ();
    my $stacksize           = $stem[2] - $stem[0] + 1;
    for ( my $x = 0 ; $x < $stacksize ; $x++ ) {

      #    for (my $i = $stem[0]; $i <= $stem[2] ; $i++){
      #      for(my $j = $stem[1]; $j >= $stem[3]; $j--){
      my $i = $stem[0] + $x;
      my $j = $stem[1] - $x;
      @p_bp = @c_bp;
      @c_bp = ( $i, $j );
      if ( defined $bps{"$i:$j"} ) {

        # add indices to stem relation object
        push( @stem_relation_idxs, $i );
        push( @stem_relation_idxs, $j );

        # add edge to graph edges
        push( @edg, "e $i $j s" );

        # add stack if option given
        if ($stacks) {
          if (@p_bp) {
            push( @stackgraph,          "v $curr_gi P" );
            push( @stackgraph,          "e $curr_gi $p_bp[0] p" );
            push( @stackgraph,          "e $curr_gi $p_bp[1] p" );
            push( @stackgraph,          "e $c_bp[0] $curr_gi p" );
            push( @stackgraph,          "e $c_bp[1] $curr_gi p" );
            push( @stack_relation_idxs, $curr_gi );

            # update graph index
            ++$curr_gi;
          }
        }
      } else {
        die "ERROR: There is no such base-pair in the base-pair hash:" . " ($i, $j)\n";
      }
    }

    #  	  }
    #    }
    # add stem to abstract graph
    # (graph indices become a bit jumbled, but saves time)
    push( @abstrgraph, "v $curr_gi ${ABSTRUCT}#S" );

    # stem object now saves the curr_gi index for that stem
    $stems{$sk} = $curr_gi;
    ++$curr_gi;

    # add relations
    push( @abstrgraph, "v $curr_gi ^R" );
    push( @abstrgraph, "e " . ( $curr_gi - 1 ) . " $curr_gi ^r" );
    foreach my $r (@stem_relation_idxs) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    if ($stacks) {
      foreach my $r (@stack_relation_idxs) {
        push( @abstrgraph, "e $curr_gi $r \@r" );
      }
    }
    ++$curr_gi;
  }

  #-------------------> add HLs
  foreach my $hl (@hls) {

    # check size
    die "ERROR: the HL object is the incorrect size\n" unless ( @{$hl} == 2 );

    # add vertex
    push( @abstrgraph, "v $curr_gi ${ABSTRUCT}#HL" );

    # add edge from stem to hl
    if ( defined( $stems{ $bps{"$hl->[0]:$hl->[1]"} } ) ) {
      push( @abstrgraph, "e " . $stems{ $bps{"$hl->[0]:$hl->[1]"} } . " $curr_gi n" );
    } else {
      die "ERROR: stem not defined in base-pair hash for $hl->[0]:$hl->[1]\n";
    }
    ++$curr_gi;

    # add relations
    push( @abstrgraph, "v $curr_gi ^R" );
    push( @abstrgraph, "e " . ( $curr_gi - 1 ) . " $curr_gi ^r" );

    # add all nodes between i to j, inclusively
    for ( my $r = $hl->[0] ; $r <= $hl->[1] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    ++$curr_gi;
  }

  #------------------> add BLs
  foreach my $bl (@bls) {

    # check size
    die "ERROR: the BL object is the incorrect size\n" unless ( scalar( @{$bl} ) == 4 );

    # add vertex
    push( @abstrgraph, "v $curr_gi ${ABSTRUCT}#BL" );

    # add edges from bl to both adjacent stems
    if ( defined( $stems{ $bps{"$bl->[0]:$bl->[1]"} } )
      && defined( $stems{ $bps{"$bl->[2]:$bl->[3]"} } ) ) {

      # add outer stem: S->BL
      push( @abstrgraph, "e " . $stems{ $bps{"$bl->[0]:$bl->[1]"} } . " $curr_gi n" );

      # add inner stem: BL->S
      push( @abstrgraph, "e $curr_gi " . $stems{ $bps{"$bl->[2]:$bl->[3]"} } . " n" );
    } else {
      die "ERROR: stems not defined for BL\n";
    }
    ++$curr_gi;

    # add relations
    push( @abstrgraph, "v $curr_gi ^R" );
    push( @abstrgraph, "e " . ( $curr_gi - 1 ) . " $curr_gi ^r" );

    # add all nodes between i to k and l to j, inclusively
    for ( my $r = $bl->[0] ; $r <= $bl->[2] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    for ( my $r = $bl->[3] ; $r <= $bl->[1] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    ++$curr_gi;
  }

  #------------------> add ILs
  foreach my $il (@ils) {

    # check size
    die "ERROR: the IL object is the incorrect size\n" unless ( scalar( @{$il} ) == 4 );

    # add vertex
    push( @abstrgraph, "v $curr_gi ${ABSTRUCT}#IL" );

    # add edges from il to both adjacent stems
    if ( defined( $stems{ $bps{"$il->[0]:$il->[1]"} } )
      && defined( $stems{ $bps{"$il->[2]:$il->[3]"} } ) ) {

      # add outer stem: S->IL
      push( @abstrgraph, "e " . $stems{ $bps{"$il->[0]:$il->[1]"} } . " $curr_gi n" );

      # add inner stem: IL->S
      push( @abstrgraph, "e $curr_gi " . $stems{ $bps{"$il->[2]:$il->[3]"} } . " n" );
    } else {
      die "ERROR: stems not defined for IL\n";
    }
    ++$curr_gi;

    # add relations
    push( @abstrgraph, "v $curr_gi ^R" );
    push( @abstrgraph, "e " . ( $curr_gi - 1 ) . " $curr_gi ^r" );

    # add all nodes between i to k and l to j, inclusively
    for ( my $r = $il->[0] ; $r <= $il->[2] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    for ( my $r = $il->[3] ; $r <= $il->[1] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    ++$curr_gi;
  }

  #------------------> add ELs
  foreach my $el (@els) {

    #check size
    die "ERROR: the EL object is the incorrect size\n" unless ( scalar( @{$el} ) == 4 );

    # add vertex
    push( @abstrgraph, "v $curr_gi ${ABSTRUCT}#EL" );

    # add edges from el to both adjacent stems, if available
    unless ( $el->[0] == $el->[1] ) {
      if ( defined( $stems{ $bps{"$el->[0]:$el->[1]"} } ) ) {
        push( @abstrgraph, "e " . $stems{ $bps{"$el->[0]:$el->[1]"} } . " $curr_gi n" );
      } else {
        die "ERROR: stems not defined for left BP of EL\n";
      }
    }
    unless ( $el->[2] == $el->[3] ) {
      if ( defined( $stems{ $bps{"$el->[2]:$el->[3]"} } ) ) {
        push( @abstrgraph, "e $curr_gi " . $stems{ $bps{"$el->[2]:$el->[3]"} } . " n" );
      } else {
        die "ERROR: stems not defined for right BP of EL\n";
      }
    }
    ++$curr_gi;

    # add relations
    push( @abstrgraph, "v $curr_gi ^R" );
    push( @abstrgraph, "e " . ( $curr_gi - 1 ) . " $curr_gi ^r" );

    # @r relation edges
    for ( my $r = $el->[1] ; $r <= $el->[2] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }

    # if BP, not EOS
    unless ( $el->[0] == $el->[1] ) {
      push( @abstrgraph, "e $curr_gi $el->[0] \@r" );
    }

    # if BP, not EOS
    unless ( $el->[2] == $el->[3] ) {
      push( @abstrgraph, "e $curr_gi $el->[3] \@r" );
    }
    ++$curr_gi;
  }

  #------------------> add MLs
  foreach my $ml (@mls) {

    #check size
    die "ERROR: the ML object is not big enough!\n" unless ( scalar( @{$ml} ) >= 3 );

    # add vertex
    push( @abstrgraph, "v $curr_gi ${ABSTRUCT}#ML" );

    # add edge from outer stem to first BP
    my $closing_bp = shift( @{$ml} );
    if ( defined( $stems{ $bps{"$closing_bp->[0]:$closing_bp->[1]"} } ) ) {
      push( @abstrgraph, "e " . $stems{ $bps{"$closing_bp->[0]:$closing_bp->[1]"} } . " $curr_gi n" );
    } else {
      die "ERROR: ML is not defined for BP $closing_bp->[0]:$closing_bp->[1]\n";
    }

    # add remaining stems
    foreach my $bp ( @{$ml} ) {
      if ( defined( $stems{ $bps{"$bp->[0]:$bp->[1]"} } ) ) {
        push( @abstrgraph, "e $curr_gi " . $stems{ $bps{"$bp->[0]:$bp->[1]"} } . " n" );
      } else {
        die "ERROR: ML is not defined for BP $bp->[0]:$bp->[1]\n";
      }
    }
    ++$curr_gi;

    # add relations
    push( @abstrgraph, "v $curr_gi ^R" );
    push( @abstrgraph, "e " . ( $curr_gi - 1 ) . " $curr_gi ^r" );

  # add @r relation edges between i and k, where (i,j) is the first BP and (k,l)
  # the second BP
    my $pbp_aref = $closing_bp;
    my $cbp_aref = shift( @{$ml} );
    for ( my $r = $pbp_aref->[0] ; $r <= $cbp_aref->[0] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }

    # rest of @r relation edges
    foreach my $bp ( @{$ml} ) {
      $pbp_aref = $cbp_aref;
      $cbp_aref = $bp;
      for ( my $r = $pbp_aref->[1] ; $r <= $cbp_aref->[0] ; $r++ ) {
        push( @abstrgraph, "e $curr_gi $r \@r" );
      }
    }

    # add final stretch between last BP (m,n) and first closing base-pair (i,j)
    for ( my $r = $cbp_aref->[1] ; $r <= $closing_bp->[1] ; $r++ ) {
      push( @abstrgraph, "e $curr_gi $r \@r" );
    }
    ++$curr_gi;
  }

  @graph_lines = ( @edg, @stackgraph, @abstrgraph );
  return ( \@graph_lines, $curr_gi );
}

# just to test the current shrep dot-bracket structure parsing
sub test {
  my ( $shrep, $idx, $gi, $mark_aref, $isML ) = @_;

  my $result = "";

  $result .= "$shrep\n";
  my @a = split( "", $shrep );
  for ( my $i = 0 ; $i <= @a ; $i++ ) {
    my $mark = 0;
    foreach my $j ( @{$mark_aref} ) {
      if ($isML) {
        $mark = 1 if ( $i == $j->[0] - $gi );
        $mark = 1 if ( scalar( @{$j} == 2 && $i == $j->[1] - $gi ) );
      } else {
        $mark = 1 if ( $i == $j - $gi );
      }

    }
    if ( $i == $idx - $gi ) {
      $result .= "^";
    } elsif ($mark) {
      $result .= "'";
    } else {
      $result .= " ";
    }
  }
  $result .= "\n";

  return $result;
}

###########################################################################
# Extends or closes a multi-loop. This means a ML already exists.
# We can have the case ((x)( or ((x).(
#
# Input:
# curr_shrep	The current shrep object, with shrep string at pos 0
# c_ml_aref		The current open MLs array reference
# p_bp_aref		The previous base-pair array reference
# stmbrks_aref	The stem-breaks array reference
# open_blks_aref	The open bracket blocks array reference
# idx			The current index
# #curr_gi		The current graph index
# Output:
# None, it should just modify the ML accordingly
############################################################################
sub extend_or_open_nested_ML {
  my ( $curr_shrep, $c_mls_aref, $p_bp_aref, $up_aref, $stmbrks_aref, $open_blks_aref, $idx, $curr_gi ) = @_;

  # case '(x(x)x(x)x(' extend current ML
  if ( $c_mls_aref->[-1]->[-1]->[-1] == $p_bp_aref->[0] ) {
    push( @{ $c_mls_aref->[-1]->[-1] }, $p_bp_aref->[1] );    #close previous BP
    my @a = ($idx);
    push( @{ $c_mls_aref->[-1] }, \@a );    # add new open base-pair
    print STDERR "TEST - extended ML:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, $c_mls_aref->[-1], 1 ) ) if ($i_debug);

    # case '(x(x)x((x)(' nested ML, => OPEN_ML with stem-break
  } elsif ( @$stmbrks_aref
    && $stmbrks_aref->[-1]->[1] == $p_bp_aref->[0]
    && $stmbrks_aref->[-1]->[2] == $p_bp_aref->[1] ) {

    my @c_stmbrk = @{ pop( @{$stmbrks_aref} ) };
    my @a        = ( $c_stmbrk[0], -1 );
    my @ml       = ();
    push( @ml, \@a );                       # 1st base-pair is in-waiting
    my @a2 = ( $c_stmbrk[1], $c_stmbrk[2] );
    push( @ml, \@a2 );                      #2nd BP
    my @a3 = ($idx);
    push( @ml,            \@a3 );           #first of current BP
    push( @{$c_mls_aref}, \@ml );
    print STDERR "TEST - opened nested ML with SB:\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@ml, 1 ) ) if ($i_debug);

    # case '(x(x)x(.(x)(' nested ML, => OPEN_ML with UP
  } elsif ( @$up_aref
    && $up_aref->[-1]->[0] - 1 == $open_blks_aref->[-1]->[-1]
    && $up_aref->[-1]->[1] + 1 == $p_bp_aref->[0] ) {
    my @tmp_up = @{ pop( @{$up_aref} ) };
    my @a      = ( $open_blks_aref->[-1]->[-1], -1 );
    my @ml     = ();
    push( @ml, \@a );                       # 1st BP awaits closing
    my @newbp = @{$p_bp_aref};
    push( @ml, \@newbp );                   #2nd BP
    my @a2 = ($idx);
    push( @ml,            \@a2 );           # current open BP, awaits closing
    push( @{$c_mls_aref}, \@ml );
    print STDERR "TEST - opened nested ML with prev UP\n" . ( test( $curr_shrep->[0], $idx, $curr_gi, \@ml, 1 ) ) if ($i_debug);
  } else {
    die "ERROR: in the case extend or open nested ML, " . "but the base-pairs don't fit for either!\n" . ( test( $curr_shrep->[0], $idx, $curr_gi ) ) if ($i_debug);
  }
}

###########################################################################
# Closes a stem-object for the abstract graph. If necessary, it also
# creates a stem-break because the current open block is still open,
# which can only mean that there is a stem-break between two consecutive
# opening brackets.
#
# Input:
# c_stem_aref : current stem array reference
# stems_href : all stems hash reference
# p_bp_aref : previous base-pair array reference
# c_cp_aref : current base-pair array reference
# p_open_aref : previously opened open-bracket block
# p_stmbrk_aref : the previous stem-break
# open_blks_aref : all old open blocks array reference
# bps_href : all base-pairs hash reference
#
# Output:
# None, it should just modify the input variables accordingly
############################################################################
sub close_stem {
  my ( $c_stem_aref, $stems_href, $p_bp_aref, $c_bp_aref, $p_open_aref, $stmbrks_aref, $open_blks_aref, $bps_href ) = @_;

  my @newstem = @$c_stem_aref;
  my $stemkey = "$newstem[0]:$newstem[1],$newstem[2]:$newstem[3]";
  $stems_href->{$stemkey} = \@newstem;

  # reset and clean running variables
  @$c_stem_aref = ();
  @$p_bp_aref   = @$c_bp_aref;
  @$c_bp_aref   = ();

  # add stem to base-pair hash
  $bps_href->{"$newstem[0]:$newstem[1]"} = $stemkey;
  $bps_href->{"$newstem[2]:$newstem[3]"} = $stemkey;

  # there are still bases open in the previous open block, '[[]x',
  # where [] is a closed stem and [ is an open block,
  # must create STEM-BREAK
  # this cannot be the case ((x)), because there p_open is empty
  if (@$p_open_aref) {

    # check: the next open base is adjacent to the previous BP
    if ( $p_open_aref->[-1] == $p_bp_aref->[0] - 1 ) {
      my @stmbrk = ( $p_open_aref->[-1], @$p_bp_aref );
      push( @$stmbrks_aref, \@stmbrk );

      # push the previous open block back onto stack
      my @newopen = @$p_open_aref;
      push( @$open_blks_aref, \@newopen );
      @$p_open_aref = ();
    } else {

      # error, this should not be possible
      die("ERROR in close_stem(), [.[]x should not occur here\n");
    }
  }
}

###########################################################################
# Here the information about the structure is added to the graph. In a
# first step each base-pair is added by adding edges between the vertices
# of the respective nucleotides. Then, if two base-pairs are stacked,
# we add an extra vertex and connect the four nucleotide vertices that
# are involved with extra edges (symbols P used for vertex and p for edges).
# While adding new vertices we keep track of the current graph index.
#
# Input:
# curr_shrep : the current shrep as dot-bracket structure
# curr_gi : the current graph index
# win_size : the window size
# stacks : the input parameter -stacks (whether to add stacking information)
#
# Output:
# The structure graph information line-by-line as an array reference and
# the current graph index.
############################################################################
sub getStructGraph {
  my ( $curr_shrep, $curr_gi, $win_size, $stacks ) = @_;

  my @struct = split( "", $curr_shrep->[0] );

  my @edg    = ();
  my @starts = ();
  my @pairs  = ();

  foreach my $idx ( 0 .. @struct - 1 ) {

    push( @starts, $idx ) if ( $struct[$idx] eq "(" );

    if ( $struct[$idx] eq ")" ) {
      my $start = $curr_gi + pop(@starts);
      my $end   = $curr_gi + $idx;
      push( @edg, "e " . $start . " " . $end . " s " );
      my @pair = ( $start, $end );
      push( @pairs, \@pair );
    }
  }

  # update current graph index
  $curr_gi = $curr_gi + $win_size;

  # initiate structure lines array
  my @stacking_info = ();

  # add stacking information to graph unless input option tells us not to
  if ($stacks) {

    my $stacked_pairs = 0;

    # get stacked base-pairs (they are ordered according to
    # position of closing nucleotide)
    for ( my $i = 1 ; $i < @pairs ; $i++ ) {

      # add stacked base-pairs (vertices+edges)
      my $curr_bp_aref = $pairs[$i];
      my $prev_bp_aref = $pairs[ $i - 1 ];

      # if the current base-pair is stacked on the previous base-pair,
      # when curr_start = prev_start - 1 AND curr_end = prev_end + 1
      if ( $curr_bp_aref->[0] == $prev_bp_aref->[0] - 1
        && $curr_bp_aref->[1] == $prev_bp_aref->[1] + 1 ) {

        # add stacking vertex P
        push( @stacking_info, "v $curr_gi P" );

        # add four edges from involved nucleotids
        push( @stacking_info, "e $curr_gi $prev_bp_aref->[0] p" );
        push( @stacking_info, "e $curr_gi $prev_bp_aref->[1] p" );
        push( @stacking_info, "e $curr_bp_aref->[0] $curr_gi p" );
        push( @stacking_info, "e $curr_bp_aref->[1] $curr_gi p" );
        ++$curr_gi;    # add one to the index, ready for next vertex

      }
    }
  }

  my @str_graphlines = ( @edg, @stacking_info );
  return ( \@str_graphlines, $curr_gi );
}

###########################################################################
# This method gives us the overall graph header for this sequence. This
# graph includes all window calculations and all shreps for all windows.
# Input:
# seq_id : the sequence ID
# seq_head : the header information for the sequence
# wins_aref : the windows sizes for folding
# h_shift : input parameter -shift
# h_e : input parameter -e
# h_c : input parameter -c
# h_t : input parameter -t
# h_u : input parameter -u
# h_r : input parameter -r
# h_M : input parameter -M
# h_crop_unpaired_ends : input parameter -cue
#
# Output:
# The header line as a string
############################################################################
sub getGraphHeader {
  my ( $seq_id, $seq_head, $wins_aref, $h_shift, $h_e, $h_c, $h_t, $h_u, $h_r, $h_M, $h_crop_unpaired_ends, $h_i, $h_sample_len, $h_q, $h_Tp, $h_seqlen, $h_noStr ) = @_;

  my $ret;

  $ret = "t # SEQID $seq_id ";
  $ret .= "$seq_head " if ( defined $seq_head );
  $ret .= "SEQLEN $h_seqlen ";

  return $ret . "\n" if ($h_noStr);

  $ret .= "MAXWINSHREPS $h_M ";
  $ret .= "CUE 1 " if ($h_crop_unpaired_ends);
  $ret .= "RNASHAPES -w " . ( join( ",", @{$wins_aref} ) ) . " ";
  $ret .= "-SHIFT%w $h_shift " if ( not $globalFolding and $h_shift );
  $ret .= "-SHIFT 1 " if ( not $globalFolding and not $h_shift );
  $ret .= "-e $h_e "                         if ( defined $h_e );
  $ret .= "-c $h_c "                         if ( defined $h_c );
  $ret .= "-t $h_t "                         if ( defined $h_t );
  $ret .= "-u 1 "                            if ($h_u);
  $ret .= "-r 1 "                            if ($h_r);
  $ret .= "-i $h_i SAMPLE_LEN $h_sample_len" if ( not $h_q and $h_i );
  $ret .= "-q 1 "                            if ($h_q);
  $ret .= "-T $h_Tp"                         if ( $h_q and defined $h_Tp );
  $ret .= "\n";

  return $ret;
}

###########################################################################
# This method gives us the window header for the subgraph that includes
# all shreps for the current window.
# Input:
# graphHead : the header line for the entire sequence graph
# win_size : the size of the current window
# win_start : the starting position of the current window
# win_end : the end position of the current window
# win_center : the centre position of the current window
# win_global_fold : (boolean) whether the actual number of nucleotides
# is smaller than the given window size (end of sequence?)
#
# Output:
# The header line as a string
############################################################################
sub getWindowHeader {
  my ( $graphHead, $win_size, $win_start, $win_end, $win_center, $win_global_fold, $win_sample, $used_t ) = @_;

  chomp $graphHead;

  $graphHead =~ s/t # //;

  my $ret = "w # $graphHead ";
  $ret .= "SHAPE_TYPE $used_t ";
  $ret .= "GLOBALFOLD $win_global_fold ";
  $ret .= "WSIZE $win_size ";
  $ret .= "WSTART $win_start ";
  $ret .= "WEND $win_end ";
  $ret .= "WCENT $win_center ";
  $ret .= "WSAMPLE $win_sample ";
  $ret .= "\n";

  return $ret;
}

###########################################################################
# This method gives us the shape header that includes the single shrep
# connected components.
# Input:
# winHead : the header for the current window
# shrep : information about the current shrep
#
# Output:
# The header line as a string
############################################################################
sub getShapeHeader {
  my ( $winHead, $shrep ) = @_;

  chomp $winHead;

  $winHead =~ s/w # //;
  $winHead =~ s/t # //;

  my $ret = "s # $winHead ";

  my @info = @{$shrep};
  $ret .= join( " ", @info[ 1 .. $#info ] );
  $ret .= "\n";

  return $ret;
}

###########################################################################
# This method gives us the sequence header that includes the unstructured
# sequence connected components.
# Input:
# winHead : the header for the current window
# shrep : information about the current shrep
#
# Output:
# The header line as a string
############################################################################
sub getSeqHeader {
  my ( $winHead, $shrep ) = @_;

  chomp $winHead;

  $winHead =~ s/w # //;
  $winHead =~ s/t # //;

  my $ret = "u # $winHead ";

  my @info = @{$shrep};
  $ret .= join( " ", @info[ 1 .. $#info ] );
  $ret .= "\n";

  return $ret;
}

1;

#############################################################################
# Programming description
#############################################################################
#	Substructure graphs for machine learning with Fabrizio
#	-------------------------------------------------------
#
#	(1) Parameters (RNAshapes parameter):
#		- Window sizes [] (-w)
#		- window shift size (-W)
#	- calculate structure probabilities (-r)
#		- energy range kcal/mol (-e) OR energy relative percentage (%) to MFE (-c)
#		- shape type 1 to 5 (-t)
#		- ignore unstable substructures (-u)
#		- max shreps
#
#
#	(2) For each sequence, generate one graph/file that consists of all windows. The general format for one graph is as follows:
#
#	t # seq_id parameters
#	v graph_index nt_type window_size window_centre abs_seq_pos
#	...
#	e m m+1 > 1 (backbone)
#	...
#	e base_i_graph_index base_j_graph_index s shrep_e shrep_p ...
#
#	For each window (subgraph) we create a subgraph (of the subgraph) for each substructure.
#	We have a running index (gi) for each subgraph. All vertex and edge indices of the subgraph add
#	the running graph index to the actual window position. For example
#
#
#	Sequence: AAACC CUUUG GG
#		  01234 56789 01
#
#	Window=10 substructure1 = (((...))). centre 5.5
#
#	v 0 A 10 5.5 0
#	v 1 A 10 5.5 1
#	v 2 A 10 5.5 2
#	v 3 C 10 5.5 3
#	v 4 C 10 5.5 4
#	v 5 C 10 5.5 5
#	v 6 U 10 5.5 6
#	v 7 U 10 5.5 7
#	v 8 U 10 5.5 8
#	v 9 G 10 5.5 9
#	e 0 1 > 1
#	e 1 2 > 1
#	e 2 3 > 1
#	e 3 4 > 1
#	e 4 5 > 1
#	e 5 6 > 1
#	e 6 7 > 1
#	e 7 8 > 1
#	e 8 9 > 1
#	e 0 8 s -15.0 0.1223
#	e 1 7 s -15.0 0.1223
#	e 2 6 s -15.0 0.1223
#
#	gi = 9+1 = 10
#
#	Window = 10 substructure2 = .(((...))) centre 6.5
#	v 10 A 10 6.5 2
#	v 11 C 10 6.5 3
#	v 12 C 10 6.5 4
#	v 13 C 10 6.5 5
#	v 14 U 10 6.5 6
#	v 15 U 10 6.5 7
#	v 16 U 10 6.5 8
#	v 17 G 10 6.5 9
#	v 18 G 10 6.5 10
#	v 19 G 10 6.5 11
#	e 10 11 > 1
#	e 11 12 > 1
#	e 12 13 > 1
#	e 13 14 > 1
#	e 14 15 > 1
#	e 15 16 > 1
#	e 16 17 > 1
#	e 17 18 > 1
#	e 18 19 > 1
#	e 11 19 s -17.0 0.156
#	e 12 18 s -17.0 0.156
#	e 13 17 s -17.0 0.156
#
#	gi = 19+1 = 20
#
#
#
#	Write one perl script to create graphs for a set of sequences, fasta2shrep_gspan.pl.
#
#	INPUT:
#		-f fasta file with all sequences to compute
#		parameters as above
#
#	OUTPUT:
#		one file per sequence that contains graph called seq_id.gspan
#
#	(1) for each window size call RNAshapes and write to a tmp file
#	(2) parse result of RNAshapes (catch RNAshapes error - sequence too long?) check for max shreps.
#	(3) convert RNAshapes result to subgraph -> write to file (readpipe) look at efficiency and errors
#	(4) repeat (1) to (3) for each sequence

# ABSTRACT STRUCTURE GRAPH
# To each shrep graph (s #) add an abstract structure graph with special
# abstract relations. For example, to add the abstract structure graph to
# the previous shrep structure with gi numbers from 10 to 19, we first
# identify the abstract shape, i.e. EL-S-HL and then create nodes (labelled
# with a given name-space followed the actual node name and separated by a
# hash) and edges for this graph as follows:
# NOTE: at the moment we add the adjacent base-pairs to the loop definitions
# but these can be removed if necessary in the future
# v 20 abstruct#EL
# v 21 abstruct#S
# v 22 abstruct#HL
# e 20 21 n
# e 21 22 n
# v 23 ^R
# e 20 23 ^r
# e 23 10 @r
# e 23 11 @r
# e 23 19 @r
# v 24 ^R
# e 21 24 ^r
# e 24 11 @r
# e 24 12 @r
# e 24 13 @r
# e 24 17 @r
# e 24 18 @r
# e 24 19 @r
# v 25 ^R
# e 22 25 ^r
# e 25 13 @r
# e 25 14 @r
# e 25 15 @r
# e 25 16 @r
# e 25 17 @r

# gi = 25+1 = 26

# ANOTATION FILE
# This is a file that labels a sequence region with a given annotation. In one
# file we can have annotations within different name-spaces, for example target
# sites predicted with different tools.
# File format is a tab-delimited file with the following columns:
# SEQID (same ID as in the fasta file header)
# i (left positon of region)
# j (right position of region, if one position i=j)
# NAMESPACE#LABEL (give each annotation type one name-space and choose a label depending on the task)

# E.g we have 2 different miRNAs and we predict the target-sites with (1) IntaRNA and (2) TargetSearch,
# then this could be IntaRNA#miR1, IntaRNA#miR2, TargetSearch#miR1, and TargetSearch#miR2.
# All labels can be used more than once for one or more sequences.
# all labels with the same namespace and the same sequence ID are grouped into
# one abstract graph, according to the order of i.
# Example
# SEQID 	i 	j 	NAMESPACE#LABEL
# s1	10 	20	IntaRNA#miR1
# s1	54	60	IntaRNA#miR2
# s1	15	25	TargetSearch#miR1
# s1	54	60 	TargetSearch#miR2
# s2 ...

# We create connected component per sequence ID per name-space,
# so that all labels with
# the same namespace are grouped together (for one sequence) as follows:
# v 26 IntaRNA#miR1
# v 27 IntaRNA#miR2
# e 26 27 n
# then create R nodes and ^r and @r relations as before (from i to j)

# If a sequence graph option is given, then these graphs (u #) have to
# be labelled with this annotation
