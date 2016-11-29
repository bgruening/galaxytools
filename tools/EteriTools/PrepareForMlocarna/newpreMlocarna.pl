#!/usr/bin/perl
use List::Util qw/ min max /;


my ( $CI, $fast_cluster, $data_fasta, $fast_cluster_sim, $map_data,$nspdk_knn_center, $myPath ) = @ARGV;



if (not defined $fast_cluster ||  not defined $data_fasta || not defined $fats_cluster_sim || not defined $map_data) {
  die "Wrong arguments: candidate clusters, fasta file, similarity scores and/or map data are missing \n";
}



#my $nspdk_knn_center = 20;
#my $ids_aref =  readSubset($fast_cluster, 2, $nspdk_knn_center);
#print $ids_aref;
#print "\n";

my $CLUSTER_DIR = "CLUSTER";
system("mkdir -p $CLUSTER_DIR");
my @fa = read_fasta_file($data_fasta);





## get real number of dense centers/clusters found by nspdk ($CONFIG{GLOBAL_num_clusters} could be higher)
my $num_clusters_curr = readpipe("wc -l $fast_cluster");
$num_clusters_curr =~ /(\d+)\s.*$/;
$num_clusters_curr = $1;

print "\n########### \n";
print  $num_clusters_curr;
print "\n###########\n";




#print "Round $CI Stage 6-8: cluster $num_clusters_curr centers with LocARNA (parallel)";

my %toDo_models = ();
map { $toDo_models{"$CI.$_"} = $_ } ( 1 .. $num_clusters_curr );

my $cluster_error = 0;

my %job_task_finished     = ();
my $trigger_new_partition = 0;


while ( keys %toDo_models ) {



    foreach my $clus_idx ( sort { $toDo_models{$a} <=> $toDo_models{$b} } keys %toDo_models ) {
#        my $model_dir = "$CLUSTER_DIR/$clus_idx.MODEL";
#        system("mkdir -p $model_dir");

  #    if ( -e "$CLUSTER_DIR/$clus_idx.cluster/cmsearch.DONE" ) {
  #        print "Round $CI cluster $clus_idx stages 6-8 already finished!\n";
  #        delete $toDo_models{$clus_idx};
  #        next;
  #    }


        $clus_idx =~ /\d+\.(\d+)/;
        my $clus_idx_ci      = $1;

        print "cluster_idx = $clus_idx \n";
        print "cluster_idx_ci = $clus_idx_ci \n";

        #my $curr_cluster_dir = "$CLUSTER_DIR/$clus_idx.cluster";

        #system("mkdir -p $CLUSTER_DIR/$clus_idx.cluster");


        my $ids_aref =  readSubset($fast_cluster, $clus_idx_ci, $nspdk_knn_center);
        writeSet( $ids_aref, "$CLUSTER_DIR/$clus_idx.center.ids" );
        writeSubsetFrags( \@fa, $ids_aref, "$CLUSTER_DIR/$clus_idx.center.frags", "SEQ" );
        writeSubsetFasta( \@fa, $ids_aref, "$CLUSTER_DIR/$clus_idx.center.fa", 1 );

        my $center_fa_file = "$CLUSTER_DIR/$clus_idx.center.fa";

        my $ids_all = readSubset( $fast_cluster, $clus_idx_ci, $nspdk_knn_center * 3 );
        writeSet( $ids_all, "$CLUSTER_DIR/$clus_idx.center.ids.ext" );
        writeSubsetFasta( \@fa, $ids_all, "$CLUSTER_DIR/$clus_idx.center.fa.ext", 1 );

        my @fa_center = read_fasta_file($center_fa_file);
        my $subset_fa = writeSubsetFasta( \@fa_center, $ids_aref, "$CLUSTER_DIR/$clus_idx.model.tree.fa", 1 );

        my $knn_sim = readSubset( $fast_cluster_sim, $clus_idx_ci );

        open( OUT, ">$CLUSTER_DIR/$clus_idx.center.sim" );
                print OUT join( " ", @{$knn_sim} ) . "\n";
        close(OUT);

        system("touch $CLUSTER_DIR/$clus_idx.pp.DONE");

        my $ids_knn = readSubset( "$CLUSTER_DIR/$clus_idx.center.ids", 1 );

        open( OUT, ">$CLUSTER_DIR/$clus_idx.names" );
        map { print OUT $_ . "\n"; } sort { $a <=> $b } @{$ids_knn};
        close(OUT);


        #my $tree_dir = "$CLUSTER_DIR/$clus_idx.TREE";
        #system("mkdir -p $tree_dir");

        sim2matrix( "$CLUSTER_DIR/$clus_idx.center.sim", $ids_knn, "$CLUSTER_DIR/$clus_idx.matrix.kernel" );


        matrix_blacklist_overlap_frags( "$CLUSTER_DIR/$clus_idx.matrix.kernel", "$CLUSTER_DIR/$clus_idx.center.ids", $map_data, "$CLUSTER_DIR/$clus_idx.matrix.kernel.ol", 0.1, 6 );


        system("cp $CLUSTER_DIR/$clus_idx.matrix.kernel.ol $CLUSTER_DIR/$clus_idx.matrix.tree");

        my $tree_file = "$CLUSTER_DIR/$clus_idx.tree";
        matrix2tree( "$CLUSTER_DIR/$clus_idx.matrix.tree", "$CLUSTER_DIR/$clus_idx.names", $CLUSTER_DIR, $tree_file );

        die "No tree file found ($tree_file)! Exit...\n\n" if ( !-e $tree_file );
        print "$center_fa_file \n";



         my $ids_ext = readSubset( "$CLUSTER_DIR/$clus_idx.center.ids.ext", 1 );
         my @fa_ext_all = read_fasta_file("$CLUSTER_DIR/$clus_idx.center.fa.ext");
         my $fa_ext_merged = mergeFrags( \@fa_ext_all );
         writeSubsetFasta( $fa_ext_merged, $fa_ext_merged->[1], "$CLUSTER_DIR/$clus_idx.cmfinder.fa", 1 );
      #   system("cp $tree_dir/mloc/results/result.aln  $model_dir/model.tree.stk ");




        delete $toDo_models{$clus_idx};

 } #end foreach

    $trigger_new_partition = 0;




}   #end while

#system("zip -r  clInfo.zip CLUSTER ");




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
	  #print $key;

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

	#print "key= $key \n";

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
  #system("cat $out_file");
  return $out_file;
}


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


sub sim2matrix {
  my $sim_file = $_[0];
  my $ids      = $_[1];    ## array ref
  my $out_file = $_[2];

  open( IN, "$sim_file" );
  my $t = <IN>;
  close(IN);
  chomp($t);
  my @ent = split( " ", $t );

  my $knn = @{$ids}; ## knn x knn matrix to build, sim file contains kernel similarities
  die "Sim matrix line has incorrect number of entries! Expected " . ( $knn * ( $knn - 1 ) / 2 ) . " found " . @ent . " Exit...\n\n"
    if ( @ent != $knn * ( $knn - 1 ) / 2 );

  ## fix order of frag ids
  my @new_order = sort { $a <=> $b } @{$ids};

  ## map frag-id to matrix col, start with col 0
  my %id2col = ();
  map { $id2col{ $new_order[ $_ - 1 ] } = $_ - 1 } 1 .. @new_order;

  ## final matrix
  my @matrix = ();
  foreach my $e (@ent) {
    my @vals = split( ":", $e );
    $matrix[ $id2col{ $vals[0] } ][ $id2col{ $vals[1] } ] = $vals[2];
    $matrix[ $id2col{ $vals[1] } ][ $id2col{ $vals[0] } ] = $vals[2];
  }

  ## diagonal
  map { $matrix[ $_ - 1 ][ $_ - 1 ] = 0 } 1 .. $knn;

#~ for($row = 0; $row < 20; $row++) {
    #~ for($col = 0; $col < 20; $col++) {
        #~ print "$matrix[$row][$col] ";
   #~ }
   #~ print "\n";
#~ }


  ## write out
  open( OUT, ">$out_file" );
  map { print OUT join( " ", @{$_} ) . "\n" } @matrix;
  close(OUT);
}


sub matrix_blacklist_overlap_frags {
  my $matrix_file    = $_[0];   ## one row per line, cols space seperated
  my $matrix_id_file = $_[1];   ## space-seperated file with all ids in one row!
  my $id_map_file = $_[2];   ## usually data.map file line: "122 SEQ23#20#100#+"
  my $matrix_out_file = $_[3];    ## filename for output matrix
  my $black_overlap = $_[4]; ## minimal required overlap of two frags to become blacklisted
  my $prec = $_[5];          ## precision for all values in matrix

  ## map id to frag
  my $frags   = read_fragments($id_map_file);
  my %id2frag = ();
  map { $id2frag{ $_->{VALUE} } = $_ } @{$frags};

  ## matrix ids
  open( IDS, $matrix_id_file );
  my @ids_matrix = <IDS>;
  close(IDS);
  chomp( $ids_matrix[0] );
  @ids_matrix = split( " ", $ids_matrix[0] );
  @ids_matrix = sort { $a <=> $b } @ids_matrix;

  ## matrix frags
  my $matrix_frags = [];
  map { push( @{$matrix_frags}, $id2frag{$_} ) } @ids_matrix;

  ## matrix frag ols
  @{$matrix_frags} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$matrix_frags};
  my $matrix_ols = fragment_overlap( $matrix_frags, $matrix_frags, $black_overlap, 1 );

  ## store all overlaps in hash
  my %idOls = ();
  foreach my $ol ( @{$matrix_ols} ) {
    $idOls{ $matrix_frags->[ $ol->[0] ]->{VALUE} . "#" . $matrix_frags->[ $ol->[1] ]->{VALUE} } = $ol->[2];
  }

  ## get min of matrix
  open( MAT, $matrix_file );
  my $matrix_min = 10**10;
  while ( my $line = <MAT> ) {
    chomp $line;
    my @ent = split( " ", $line );

    map { $matrix_min = $_ if ( $_ < $matrix_min && $_ != 0 ) } @ent;
  }
  close(MAT);
  $matrix_min = 0 if ( $matrix_min > 0 && $matrix_min < 1 );
  #print "matrix min $matrix_min\n";

  my $matrix_blacklist_value = $matrix_min - abs( $matrix_min / 2 );

  ## create new matrix
  my @newMAT = ();
  open( MAT, $matrix_file );

  my $row = 0;
  while ( my $line = <MAT> ) {

    chomp $line;
    my @ent = split( " ", $line );

    foreach my $col ( 0 .. $#ent ) {

      ## check if for current (row,col) exists overlap
      if ( exists $idOls{ $ids_matrix[$row] . "#" . $ids_matrix[$col] } || exists $idOls{ $ids_matrix[$col] . "#" . $ids_matrix[$row] } ) {
        $ent[$col] = $matrix_blacklist_value;
      } else {
        $ent[$col] = sprintf( "%1." . $prec . "f", $ent[$col] );
      }

    }
    push( @newMAT, join( " ", @ent ) );
    $row++;
  }
  close(MAT);

  ## write new matrix
  open( OUT, ">$matrix_out_file" );
  map { print OUT $_ . "\n" } @newMAT;
  close(OUT);
  makeCol($matrix_out_file);

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

sub matrix2tree {
  my $matrix_file  = $_[0];
  my $names_file   = $_[1];
  my $tree_dir     = $_[2];
  my $tree_outfile = $_[3];

  my $temp_file;

  #my $path = "/home/eteri/GalaxyProject/galaxy/tools/PrepareForMlocarna/";
  #mkdir($tree_dir);
  system("cat $matrix_file | awk '{for(i=1;i<NR;i++){print NR,i,\$(i)}}' > $tree_dir/tree.score-list");

  system("pwd");

  system("perl $myPath/rnaclustScores2Dist.pl --quantile 1.0 < $tree_dir/tree.score-list > $tree_dir/tree.dist-list") == 0 or die " .1. command was unable to run to completion:\n\n";

  system("$myPath/./pgma $names_file $tree_dir/tree.dist-list > $tree_outfile") == 0 or die " .2. command was unable to run to completion:\n\n";

  #system(" cd $path && ./pgma $names_file $tree_dir/tree.dist-list > $tree_outfile && cd -") == 0 or die "command was unable to run to completion:\n\n";
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
    $frags{$id}               = str2frag($1);
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

      my $ols = fragment_overlap( [ $fr[0] ], \@fr, 0.1, 0 );
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
