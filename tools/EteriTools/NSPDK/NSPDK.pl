#!/usr/bin/perl
use List::Util qw/ min max /;
use POSIX qw(ceil floor);
use Array::Utils qw(:all);



my $CI = 1; ## iteration num
my $SVECTOR_DIR = "SVECTOR";
my $nspdk_mi = 0;



my ($data_fasta, $gspan, $rad, $dist, $noCache, $ensf,  $oc, $usn, $nspdk_knn_center, $nhf) = @ARGV;

#print $nspdk_knn_center;
my $group_gspan = "group.gspan";
#my $group_gspan = $gspan;
#print "gspan = $gspan \n";
#print system("file $gspan") ;
#print "\n";
#print system("file $data_fasta") ;

system("bzcat $gspan > $group_gspan");


my @fa  = read_fasta_file($data_fasta);
my $num_seqs = @{ $fa[1] };
print "\nNumber of sequences in FASTA/data.fasta: " . $num_seqs . "\n";

die "Fasta file $FASTA_DIR/data.fasta contains only $num_seqs sequences! Exit...\n\n" if ( $num_seqs <= 2 );



system("NSPDK -R $rad -D $dist -gt DIRECTED -fg $group_gspan -of");
system("mkdir -p $SVECTOR_DIR");
system("cat $group_gspan.feature_bin > $SVECTOR_DIR/data.svector.$CI");
system("rm $group_gspan.feature_bin");
system("rm $group_gspan.feature");


my $clusters_last_round = [];
if ( $CI > 1 ) {
    $clusters_last_round = foundClusters( $CI - 1 );
    print $clusters_last_round;
    my $min_diff = 0;

    my $center_last_round = 0;

    if ( -e "$SVECTOR_DIR/data.svector." . ( $CI - 1 ) . ".fast_cluster" ) {
          my @t = readpipe( "cat $SVECTOR_DIR/data.svector." . ( $CI - 1 ) . ".fast_cluster" );
          chomp(@t);
          $center_last_round = @t;
          print " aaa";
        }


        ### this part is not done yet, but because we do now only 1 round we can skip ti for now
}

## we always have blacklist paramter for NSPDK stage 5 -> need always a file
system("touch $SVECTOR_DIR/data.svector.blacklist.$CI") if ( !-e "$SVECTOR_DIR/data.svector.blacklist.$CI" );
## todo: !!! tempBlacklist is not working, needs refactor but currently we can live without it
system("touch $SVECTOR_DIR/blacklist.no_model");

####here should also be part for more rounds. check later

my $blacklist_curr_size = 0;
my @bl_num = readpipe("cat $SVECTOR_DIR/data.svector.blacklist.$CI");
chomp(@bl_num);
$blacklist_curr_size = unique(@bl_num);
print "\nFINAL blacklist size for round $CI: $blacklist_curr_size\n";

## $in_stage_end=10 -> special mode: only use NSPDK and predict candidate cluster,
## but do nothing else, stop each iteration after stage 5
if ( -e "$SVECTOR_DIR/data.svector.blacklist.$CI.special" ) {
    system("cp $SVECTOR_DIR/data.svector.blacklist.$CI.special $SVECTOR_DIR/data.svector.blacklist.$CI");
    $in_stage_end = 10;
    $blacklist_curr_size = `cat $SVECTOR_DIR/data.svector.blacklist.$CI.special | wc -l`;
}

my $num_seqs_left = $num_seqs - $blacklist_curr_size;
if ( $num_seqs_left < 2 ) {
    ## do not call NSPDK again
    print "\n...only $num_seqs_left sequences left for clustering.\n";
    print " Skip all following iterations!\n\n";

    #~ system("echo  \`date\` >> $EVAL_DIR/times/time.round.$CI"); ## only for time measurment
    #~ system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" ); ## only for time measurment
    last;
}



$OPTS_nspdk_centers = "-ensf $ensf $oc $usn";

 if ( !-e "$SVECTOR_DIR/data.svector.fast_cluster.$CI.DONE" ) {

    my $job_name = "stage 5 (fastCluster.NSPDK)";
    my $job_uuid = "stage5-$CI";

    ## fix nspdk binsize to 1000 instances to be independent of dataset size
    $OPTS_nspdk_centers =~ s/-msb\s+\S+//;
    my $nspdk_max_binsize = 1000;
    if ( $num_seqs_left > $nspdk_max_binsize ) {
      my $msb = sprintf( "%.4f", ( $nspdk_max_binsize / $num_seqs_left ) );
      $msb = 0.0001 if ( $msb < 0.0001 );
      $OPTS_nspdk_centers .= " -msb $msb";
      print "use msb $msb (num seqs left:$num_seqs_left binsize:$nspdk_max_binsize)\n";
    } else {
      $OPTS_nspdk_centers .= " -msb 1";
    }

    ## nspdk uses vector filename for output (append output type like .fast_cluster .knn .approx_knn)
    ## we need this here to use different names in case of stage_end=5,
    ## i.e. there could be multiple running nspdk instances on the same feature vector

    #system("ln -f -s $SVECTOR_DIR/data.svector $SVECTOR_DIR/data.svector.$CI");


    ### don't need this part
    my $CMD_fastClusterNSPDK = [];
    $CMD_fastClusterNSPDK->[0] = "$BIN_DIR/NSPDK";
    $CMD_fastClusterNSPDK->[1] =
      "-no-cache -rs $CI -fsb $SVECTOR_DIR/data.svector.$CI -bl $SVECTOR_DIR/data.svector.blacklist.$CI $OPTS_nspdk_centers -knn $nspdk_knn_center -ss " . $GLOBAL_num_clusters . " -nhf $nspdk_nhf -mi $nspdk_mi -fcs $nspdk_fcs ";


    ## add special options for greylist clustering
    if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {
      #my $greylist_size = `cat $FASTA_DIR/$DATA_prefix.greylist | wc -l`;
      #chomp($greylist_size);
      $CMD_fastClusterNSPDK->[2] = " -gl $FASTA_DIR/$DATA_prefix.greylist -fcs 1 -otknn "; # -knn ".($nspdk_knn_center+$greylist_size)." ";
    }

    ## add special debug options for NSPDK
    #if ($in_debug) {
    #  $CMD_fastClusterNSPDK->[3] = " -oaknn ";
    #}

    ## get size of feature vector and estimate memory requirement
    my $size = -s "$SVECTOR_DIR/data.svector";
    $size = ceil( $size / 1000000000 ) + ceil( $nspdk_nhf * 0.005 );
    $size = $size * ( $num_seqs_left / $num_seqs ) + 1;
    $size = sprintf( "%.1f", $size );
    print "used qsub h_vmem=$size" . "G\n";
    print "use threads $NUM_THREADS use SGE-PE threads $SGE_PE_THREADS\n";
    my $qsub_opts = " -l h_vmem=$size" . "G ";
    $qsub_opts .= " -v OMP_STACKSIZE=20M  -v OMP_NUM_THREADS=$SGE_PE_THREADS ";
    $qsub_opts .= " -pe \"$in_SGE_PE_name\" 1-$SGE_PE_THREADS " if ($SGE_PE_THREADS>1); ## use '-hard -l c_op2356=1' to have always same cpu model (FREIBURG only)

    ## we wait until NSPDK is finished, do not wait in special blacklist case == 10 (only for internal debug)
    #~ my $sge_wait = 1;
    #~ $sge_wait = 0 if ( $in_stage_end == 10 );

    #~ my $sge_status =
      #~ job_call( $job_name, "$BIN_DIR/fastCluster.NSPDK.sge", $CMD_fastClusterNSPDK, 1, $SGE_ERR_DIR, $in_USE_SGE, "$SVECTOR_DIR/sge_log_stage5.$CI", "$EVAL_DIR/times/time.stage.5.$CI", $sge_wait, $qsub_opts, 1, $job_uuid );

    system("NSPDK $noCache -rs $CI -fsb $SVECTOR_DIR/data.svector.$CI -bl $SVECTOR_DIR/data.svector.blacklist.$CI $OPTS_nspdk_centers -knn $nspdk_knn_center -ss " . 100 . " -nhf $nhf -mi $nspdk_mi -fcs 1");


      ## prefilter centers for interesting ones,
      ## i.e. centers which contain non greylist ids
      if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {
        system("cp $SVECTOR_DIR/data.svector.$CI.fast_cluster $SVECTOR_DIR/data.svector.$CI.fast_cluster_orig");
        ## filterGreylistCenters("$FASTA_DIR/$DATA_prefix.greylist","$SVECTOR_DIR/data.svector.$CI.fast_cluster",$nspdk_knn_center);
        filterGreylistCenters( "$FASTA_DIR/$DATA_prefix.greylist", "$SVECTOR_DIR/data.svector.$CI.approx_knn", $nspdk_knn_center );
        system("cp $SVECTOR_DIR/data.svector.$CI.approx_knn.gl_filter $SVECTOR_DIR/data.svector.$CI.fast_cluster");
      }
      #system("rm $SVECTOR_DIR/data.svector.$CI");
      system("touch $SVECTOR_DIR/data.svector.fast_cluster.$CI.DONE");

      if ($evaluate) {
        system_call(
          "$BIN_DIR/clusterQuality.sh $FASTA_DIR/class.hash $SVECTOR_DIR/data.svector.$CI.fast_cluster $SVECTOR_DIR/data.svector.$CI.fast_cluster $nspdk_knn_center $EVAL_DIR/svector $FASTA_DIR/class.size $FASTA_DIR/$DATA_prefix.map > $EVAL_DIR/svector/$CI.centers_qual",
          $in_verbose
        );
        system_call("mv $EVAL_DIR/svector/cluster.counts $EVAL_DIR/svector/$CI.cluster.counts");
        system_call("mv $EVAL_DIR/svector/cluster.class.all $EVAL_DIR/svector/$CI.cluster.class.all");
        system_call("mv $EVAL_DIR/svector/cluster.class $EVAL_DIR/svector/$CI.cluster.class");
        system_call("mv $EVAL_DIR/svector/cluster.class.frags $EVAL_DIR/svector/$CI.cluster.class.frags");
        system_call("mv $EVAL_DIR/svector/class.size $EVAL_DIR/svector/$CI.class.size");

        GraphClust::evalSVECTOR( "$EVAL_DIR/svector", $CI, "$EVAL_DIR/stage5.svector_qual.final" );
      }



  } else {
    print "Round $CI fast_cluster already done\n";
  }

#system("rm 1.group.gspan");
#system("rm 1.group.gspan.feature");

sub foundClusters {
  my $round = $_[0];

  my $part = read_partition("$EVAL_DIR/partitions/$round.soft");

  my %clusters = ();

  map { $clusters{$1} = 1 if ( $_->[3] =~ /^$round\.(\d+)$/ ) } @{$part};

  my @res = sort { $a <=> $b } keys %clusters;

  return \@res;
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
