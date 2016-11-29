#!/usr/bin/perl
use List::Util qw/ min max /;
use POSIX qw(ceil floor);
use Array::Utils qw(:all);



#my $CI = 1; ## iteration num
my $SVECTOR_DIR = "SVECTOR";
my $DATA_prefix = "data";
my $GLOBAL_hit_blacklist_overlap = 0.2; #add as param
#my $GLOBAL_num_clusters = 100;

#my $data_fasta = "FASTA/data.fasta";
$final_partition_soft="";
my $bl_list ="";
my $CI=1;
my ( $data_fasta, $data_names, $noCache, $ensf,  $oc, $usn, $nspdk_knn_center, $nhf, $nspdk_nhf_max, $nspdk_nhf_step, $GLOBAL_num_clusters, $path, $CI, $bl_list, $final_partition_soft, $fast_cluster_last_round, $GLOBAL_hit_blacklist_overlap) = @ARGV;

if ($CI == ""){
  $CI=1;
}

my $data_svector = "$SVECTOR_DIR/data.svector";


my $nspdk_mi_max  = int($nspdk_knn_center/2);
my $nspdk_mi_step = max(1, int($nspdk_mi_max/5)) ;
my $nspdk_mi      = 0;

print "data_svector = $data_svector \n";
print "data_fasta = $data_fasta \n";
print "data_names = $data_names \n";
print "noCache = $noCache \n";
print "ensf = $ensf \n";
print "oc = $oc \n";
print "usn = $usn \n";
print "nspdk_knn_center = $nspdk_knn_center \n";
print "nhf = $nhf \n";
print "nspdk_nhf_max = $nspdk_nhf_max \n";
print "nspdk_nhf_step = $nspdk_nhf_step \n";
print "CI = $CI \n";
print "final_partition_soft = $final_partition_soft \n";
#print "round_soft = $round_soft \n";
print "fast_cluster_last_round = $fast_cluster_last_round \n";
print "blacklist = $bl_list \n";
print "GLOBAL_hit_blacklist_overlap = $GLOBAL_hit_blacklist_overlap \n";
print "GLOBAL_num_clusters = $GLOBAL_num_clusters \n";
print "path = $path \n";




my @fa  = read_fasta_file($data_fasta);
my $num_seqs = @{ $fa[1] };

#system("mkdir -p $SVECTOR_DIR");

my $clusters_last_round = [];

if ( $CI > 1 ) {



  ## get list of center-idx from last round (soft partition) which lead to cluster (> results_min_cluster_size)
  $clusters_last_round = foundClusters( $CI - 1 );
  my $min_diff = 0;

  #SUBSECTION("Parameters for round $CI");

  my $center_last_round = 0;
  ## fast_cluster file not exists if $in_stage_end = 5 (also used with special blacklist)

  if ( -e $fast_cluster_last_round ) {
    my @t = readpipe( "cat $fast_cluster_last_round" );
    chomp(@t);
    $center_last_round = @t;
    print "center last round = $center_last_round \n";
  }

  print "\nModels found in last round: " . @{$clusters_last_round} . "\n";
  if ( @{$clusters_last_round} ) {
    print "\n" . join( "\n", map { ( $CI - 1 ) . ".$_.cluster/MODEL" } @{$clusters_last_round} ) . "\n";
  }

  my $tr1 = 0;

  if ( $center_last_round <= $GLOBAL_num_clusters / 2 && $nspdk_mi + $nspdk_mi_step <= $nspdk_mi_max && $nspdk_mi_step > 0 ) {
    $nspdk_mi += $nspdk_mi_step;
    print "\nToo few dense centers found in last round!\n";
    print "Set new overlap of dense regions: $nspdk_mi\n";
  } elsif ( !@{$clusters_last_round} ) { $tr1 += 1 }

  if ( @{$clusters_last_round} <= $center_last_round * (3/5) && $nhf + $nspdk_nhf_step <= $nspdk_nhf_max && $nspdk_nhf_step > 0 ) {
    $nhf += $nspdk_nhf_step;
    print "\nOnly few clusters found in last round!\n";
    print "Set new number of hash functions: $nhf\n";
  } elsif (!@{$clusters_last_round}) { $tr1 += 1 }

  if ( $tr1 >= 2 ) {
    print "No cluster found in last round!\n";
    print "\nDense region overlap and number of hash functions is at maximum! (mi=$nspdk_mi, nhf=$nhf )\n";
    print "No new iteration is started!\n\n";
    #SECTION("Trigger set to stop terative GraphClust now!");
  #  system("echo  \`date\` >> $EVAL_DIR/times/time.round.$CI"); ## only for time measurment
  #  system( "echo " . time . " >> $EVAL_DIR/times/time.round.$CI" ); ## only for time measurment
    last;
  }

} ## end if ( $CI > 1 )


      ## we always have blacklist paramter for NSPDK stage 5 -> need always a file
      system("touch $SVECTOR_DIR/data.svector.blacklist.$CI") if ( !-e "$SVECTOR_DIR/data.svector.blacklist.$CI" );
      ## todo: !!! tempBlacklist is not working, needs refactor but currently we can live without it
      system("touch $SVECTOR_DIR/blacklist.no_model");
      $bl_list = "$SVECTOR_DIR/data.svector.blacklist.$CI";


####here should also be part for more rounds. check later

if ( $CI > 1 && !-e "$SVECTOR_DIR/data.svector.blacklist.$CI.special" ) {
  ## make blacklist for this round from all Infernal hits so far
  # collect all cmsearch hits from last round as blacklist
  # hits are based on hard merged partition file, i.e. on last glob_results.pl call
  makeBlacklist( $CI, "$SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits", $GLOBAL_hit_blacklist_overlap );
}

if ( $CI > 1 ) {
  #system( "cat $SVECTOR_DIR/blacklist.no_model $SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits $FASTA_DIR/$DATA_prefix.no_match > $SVECTOR_DIR/data.svector.blacklist.$CI" );
  #system( "cat $SVECTOR_DIR/blacklist.no_model $SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits > $SVECTOR_DIR/data.svector.blacklist.$CI" );
  system( "cat  $SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits > $bl_list" );
} else {
  #system("cat $SVECTOR_DIR/blacklist.no_model $FASTA_DIR/$DATA_prefix.no_match > $SVECTOR_DIR/data.svector.blacklist.$CI");
  system("cat $SVECTOR_DIR/blacklist.no_model  > $bl_list");
}


#################################################################


my $blacklist_curr_size = 0;
my @bl_num = readpipe("cat $bl_list");
chomp(@bl_num);
$blacklist_curr_size = unique(@bl_num);
print "\nFINAL blacklist size for round $CI: $blacklist_curr_size\n";

# $in_stage_end=10 -> special mode: only use NSPDK and predict candidate cluster,
## but do nothing else, stop each iteration after stage 5
if ( -e "$SVECTOR_DIR/data.svector.blacklist.$CI.special" ) {
    system("cp $SVECTOR_DIR/data.svector.blacklist.$CI.special $SVECTOR_DIR/data.svector.blacklist.$CI");
    $in_stage_end = 10;
    $blacklist_curr_size = `cat $SVECTOR_DIR/data.svector.blacklist.$CI.special | wc -l`;
 }

## do nothing for this round if we have to few sequences left
my $num_seqs_left = $num_seqs - $blacklist_curr_size;
print "num of seqs = $num_seqs \n";
print "blacklist_curr_size = $blacklist_curr_size \n";
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

   print "stex piti mtni hastat \n";
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

    #print "options =  $OPTS_nspdk_centers \n";
    ## nspdk uses vector filename for output (append output type like .fast_cluster .knn .approx_knn)
    ## we need this here to use different names in case of stage_end=5,
    ## i.e. there could be multiple running nspdk instances on the same feature vector


##!!!!!!!!!!!!!!!!!!!!  working!!!!!!!!!!!!!!###############

    system("cd $SVECTOR_DIR && ln -f -s  data.svector data.svector.$CI");
    #link("$data_svector", "$data_svector.$CI");
    #symlink("$data_svector","$data_svector.$CI") || die "cannot symlink ";


      ## add special options for greylist clustering
      if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {

        $CMD_fastClusterNSPDK->[2] = " -gl $FASTA_DIR/$DATA_prefix.greylist -fcs 1 -otknn "; # -knn ".($nspdk_knn_center+$greylist_size)." ";
      }
####es ifi mase petqa erevi avelacnel es sistemi mej

  system("$path./NSPDK $noCache -rs $CI -fsb $data_svector.$CI -bl $bl_list $OPTS_nspdk_centers -knn $nspdk_knn_center -ss $GLOBAL_num_clusters -nhf $nhf -mi $nspdk_mi -fcs 1");
  print "NKSPDKIC heto! \n";
    #system("NSPDK $noCache -rs $CI -fsb $data_svector.$CI -bl $SVECTOR_DIR/data.svector.blacklist.$CI $OPTS_nspdk_centers -knn $nspdk_knn_center -ss $GLOBAL_num_clusters -nhf $nhf -mi $nspdk_mi -fcs 1");
    #print "data_fast = $data_svector.$CI.fast_cluster \n";

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
      system("cd $SVECTOR_DIR && rm data.svector.$CI ");

  } else {
    print "Round $CI fast_cluster already done\n";
  }

#system("rm 1.group.gspan");
#system("rm 1.group.gspan.feature");

sub foundClusters {
  my $round = $_[0];

  #my $part = read_partition("$RESULTS/partitions/$round.soft");
  #my $part = read_partition("partitions/1.soft");
  my $part = read_partition("$final_partition_soft");

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

sub makeBlacklist {
  my $curr_ci   = $_[0];
  my $bl_name   = $_[1];
  my $bl_min_ol = $_[2];

  my @blacklist = ();

  ## not necessary !?
#  system_call( "$BIN_DIR/glob_results.pl --root $in_ROOTDIR --all --summary-only --last-ci " . ( $curr_ci - 1 ), $in_verbose );

my $fragsDATA = read_fragments("$data_names");
my $finalPart = read_partition("$final_partition_soft");

#my $finalPart = read_partition("partitions/final_partition.soft");
  #my $fragsDATA = read_fragments("$FASTA_DIR/$DATA_prefix.names");
  #my $finalPart = read_partition("$RESULTS_DIR/partitions/final_partition.soft"); ###xml nenc grel vor arajini jamanak lini sovorakan, heto vorpes argument stana results.

  my $finalHits = [];
  map { push( @{$finalHits}, $_->[0] ) } @{$finalPart};
  my $fragsCMSEARCH = list2frags($finalHits);

  ## we blacklist for all hits fragments on reverse strand of hit as well!
  ## this is probably the only situation where we can ignore strand in
  ## GraphClust::fragment_overlap(); all other overlap checks are strand aware!
  ## $bl_min_ol can be larger as we base all overlap on shorter fragment in ol check
  my $ignore_strand = 1;
  my $overlaps = fragment_overlap( $fragsCMSEARCH, $fragsDATA, $bl_min_ol, $ignore_strand );

  foreach my $ol ( @{$overlaps} ) {
    push( @blacklist, $fragsDATA->[ $ol->[1] ]->{VALUE} );
  }

  my @uniq_bl = unique(@blacklist);

  my %gl = ();
  if ( -e "$FASTA_DIR/$DATA_prefix.greylist" ) {
    open( GL, "$FASTA_DIR/$DATA_prefix.greylist" );
    map { my $key = $_; chomp($key); $gl{$key} = 1 } <GL>;
    close(GL);
  }

  open( OUT, ">$bl_name" );
  foreach my $key ( sort { $a <=> $b } @uniq_bl ) {
    ## write all out but do not blacklist greylist ids, can be used in non-greylist mode as well
    print OUT $key . "\n" if ( !exists $gl{$key} );
  }
  close(OUT);

  print "\nblacklist size (infernal hits) for round $curr_ci: " . @uniq_bl . "\n";
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


sub filterGreylistCenters {
  my $gl_file     = $_[0];
  my $center_file = $_[1];
  my $knn_final   = $_[2];

  my %gl = ();
  open( GL, $gl_file );
  map { my $key = $_; chomp($key); $gl{$key} = 1 } <GL>;
  close(GL);

  open( FC, $center_file );
  my @centers = <FC>;
  chomp(@centers);
  close(FC);

  foreach my $idx ( 0 .. $#centers ) {
    my @tmp = split( " ", $centers[$idx] );
    $centers[$idx] = \@tmp[ 0 .. ( $knn_final - 1 ) ];
  }

  my @newCenters = ();

  foreach my $idx ( 0 .. $#centers ) {

    my %hits = ();
    map { $hits{$_} = 1 if ( !exists $gl{$_} ) } @{ $centers[$idx] };

    next if ( ( keys %hits ) == 0 );

    my @newCent = @{ $centers[$idx] };


    push( @newCenters, \@newCent );
  }

  open( OUT, ">$center_file.gl_filter" );

  foreach my $cent (@newCenters) {
    print OUT join( " ", @{$cent} ) . "\n";
  }

  close(OUT);

  #  system("cp $center_file $center_file.orig");
  #  system("cp $center_file.gl_filter $center_file");
}
