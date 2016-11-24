use List::Util qw/ min max /;


my ($center_fa_file, $tree_file, $p, $max_diff_am, $tau, $max_diff, $path) = @ARGV;


$mloc_opts = "-p $p --max-diff-am $max_diff_am --tau $tau  --max-diff $max_diff --alifold-consensus-dp";

#$num = $center_fa_file;
#$num =~ s/\D+\z//;
#print $str;
#print "\n";

print "tree file = $tree_file\n";

my $cmd_stk = "perl $path/mloc2stockholm.pl -file $tree_file";
$cmd_stk .= " -split_input yes -con_struct $tree_file.alifold";

#system("$cmd_stk");# == 0 or die " mloc2stockholm command was unable to run to completion:\n\n";

my $model_dir = ".";
my $dp_dir = ".";
#mlocarna_center( $center_fa_file, $model_dir, $dp_dir, 0 );

system("mlocarna $mloc_opts --treefile $tree_file $center_fa_file --tgtdir mloc") == 0 or die " mloc command was unable to run to completion:\n\n";




sub mlocarna_center {
  my $fasta    = $_[0];
  my $dir      = $_[1];
  my $dpDir    = $_[2];
  my $use_locP = $_[3];

#  my $OPTS_locarna_p_model = $CONFIG{OPTS_locarna_p_model};
#  my $OPTS_locarna_model   = $CONFIG{OPTS_locarna_model};
#  my $loc_path             = $CONFIG{PATH_LOCARNA};

  #  my $vrna_path            = $CONFIG{PATH_VRNA};

  my $loc_pp_dir = "$dir/input";
  system("mkdir -p $loc_pp_dir");

  my @fa = read_fasta_file($fasta);
  foreach my $key ( keys %{ $fa[0] } ) {
    system("ln -f -s $dpDir/$key $loc_pp_dir/$key") if ( -e "$dpDir/$key" && !$use_locP );
    print "key = $key \n";
  }

  if ($use_locP) {

    print "ifi  mej \n";

    ## check mlocarna internal plfold call opts, like prob-cutoff etc.
    ## really use plfold? prob-cutoff is not set in mlocarna and default is 0.01
    system( "mlocarna $OPTS_locarna_p_model --verbose --probabilistic --tgtdir $dir $fasta > $dir/locarnaP.out 2>&1 ", 1 );
    my $rel_cmd = "reliability-profile.pl -v -fit-once-on --structure-weight=1  --fit-penalty 0.01 --beta 200 --show-sw --out=$dir/mloc/results/result.aln.rel_plot.pdf $dir/mloc ";
    print $rel_cmd. "\n";
    my @rel = readpipe($rel_cmd);
    open( REL, ">$dir/mloc/results/result.aln.rel_signal " );
    print REL @rel[ 0 .. 2 ];
    close(REL);

  } else {

    print "elsi mej \n";

   #system("mlocarna --sequ-local true --tgtdir $dir $fasta $in_mlocarna_opts");
    system( "mlocarna $OPTS_locarna_model --verbose --skip-pp --tgtdir $dir $fasta > $dir/locarna.out 2>&1 ", 1 );
  }


  ## cleanup
#  system("rm -r -f $dir/input $dir/intermediates $dir/probs $dir/results/single_reliabilities $dir/results/ali* $dir/alifold.out");

}




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




#my $ids_ext = readSubset( "$CLUSTER_DIR/$clus_idx.cluster/center.ids.ext", 1 );
#my @fa_ext_all = read_fasta_file("$CLUSTER_DIR/$clus_idx.cluster/center.fa.ext");
#my $fa_ext_merged = mergeFrags( \@fa_ext_all );
#writeSubsetFasta( $fa_ext_merged, $fa_ext_merged->[1], "$model_dir/cmfinder.fa", 1 );
#system("cp $tree_dir/mloc/results/result.aln  $model_dir/model.tree.stk ");
