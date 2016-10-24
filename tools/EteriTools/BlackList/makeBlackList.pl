#!/usr/bin/perl
use List::Util qw/ min max /;
use POSIX qw(ceil floor);
use Array::Utils qw(:all);

#makeBlacklist( $CI, "$SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits", $GLOBAL_hit_blacklist_overlap );

my $CI =  $ARGV[0];
my  $data_names = "";
my  $final_partition_soft = "";
my $GLOBAL_hit_blacklist_overlap ="";
my $oldBlackList = "";
if( $CI > 1){
    print "ifi me \n";
    $data_names = $ARGV[1];
    $final_partition_soft = $ARGV[2];
    $oldBlackList = $ARGV[3];
    $GLOBAL_hit_blacklist_overlap = $ARGV[4];
}

print "ci = $CI \n";
print "data_names = $data_names \n";
print "final_partition_soft = $final_partition_soft \n";
print "GLOBAL_hit_blacklist_overlap = $GLOBAL_hit_blacklist_overlap \n";


#my ($CI, $data_names, $final_partition_soft, $GLOBAL_hit_blacklist_overlap )= @ARGV;

my $SVECTOR_DIR = "SVECTOR";
system("mkdir -p $SVECTOR_DIR");

system("touch $SVECTOR_DIR/data.svector.blacklist.$CI") if ( !-e "$SVECTOR_DIR/data.svector.blacklist.$CI" );
## todo: !!! tempBlacklist is not working, needs refactor but currently we can live without it
system("touch $SVECTOR_DIR/blacklist.no_model");


if ( $CI > 1 && !-e "$SVECTOR_DIR/data.svector.blacklist.$CI.special" ) {
  ## make blacklist for this round from all Infernal hits so far
  ## collect all cmsearch hits from last round as blacklist
  ## hits are based on hard merged partition file, i.e. on last glob_results.pl call
  makeBlacklist( $CI, "$SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits", $GLOBAL_hit_blacklist_overlap );
}


if ( $CI > 1 ) {
  #system( "cat $SVECTOR_DIR/blacklist.no_model $SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits $FASTA_DIR/$DATA_prefix.no_match > $SVECTOR_DIR/data.svector.blacklist.$CI" );
  system( "cat $oldBlackList $SVECTOR_DIR/round." . ( $CI - 1 ) . ".hits > $SVECTOR_DIR/data.svector.blacklist.$CI" );
} else {
  #system("cat $SVECTOR_DIR/blacklist.no_model $FASTA_DIR/$DATA_prefix.no_match > $SVECTOR_DIR/data.svector.blacklist.$CI");
  system("cat $SVECTOR_DIR/blacklist.no_model  > $SVECTOR_DIR/data.svector.blacklist.$CI");
}

my $blacklist_curr_size = 0;
my @bl_num = readpipe("cat $SVECTOR_DIR/data.svector.blacklist.$CI");
chomp(@bl_num);
$blacklist_curr_size = unique(@bl_num);
print "\nFINAL blacklist size for round $CI: $blacklist_curr_size\n";


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


  open( OUT, ">$bl_name" );
  foreach my $key ( sort { $a <=> $b } @uniq_bl ) {
    ## write all out but do not blacklist greylist ids, can be used in non-greylist mode as well
    print OUT $key . "\n" if ( !exists $gl{$key} );
  }
  close(OUT);

  print "\nblacklist size (infernal hits) for round $curr_ci: " . @uniq_bl . "\n";
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
