#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use POSIX qw(ceil floor);
use Cwd qw(abs_path getcwd);
use File::Path;
use Data::Dumper;
use Array::Utils qw(:all);
use List::Util qw/ min max /;

my $in_eval_mode;
my $in_verbose = 0;
my $in_root_dir;
my $results_top_num = 10;
$in_root_dir = "";

my ($part_type, @modTreeFiles) = @ARGV;
my $part_file;

if($part_type == 0){
   $part_file = "RESULTS/partitions/final_partition.hard.merged";
}
else{

  $part_file = "RESULTS/partitions/final_partition.soft";

}
## summary contains evaluation info for used partition
my $summary;
## part: array of [SEQ700#26#65#+,42.40,-,1.2,1.1,1,63,63]
my $part = read_partition($part_file);
my @res_clus_sort = sort { $a <=> $b } unique( map { $_->[5] } @{$part} );
my $myFasta = "FASTA/data.fasta.scan";
my @fa_scan = read_fasta_file($myFasta);
my @res_todo = ();
push( @res_todo, ( 1 .. @res_clus_sort ) );

foreach my $res_idx (@res_todo) {
  my %clus_hits = ();
  my $clus_idx  = $res_clus_sort[ $res_idx - 1 ];

  foreach my $p ( grep { $_->[5] eq $clus_idx } @{$part} ) {
    my $key = $p->[0];
    $clus_hits{$key}         = {};
    $clus_hits{$key}->{TYPE} = "CMSEARCH";
    $clus_hits{$key}->{PART} = $p;
    $clus_hits{$key}->{KEY}  = $key;
  }

  print "Warning! Used partition $part_file for cluster $clus_idx contains multiple fragments with same location!\n\n"
    if ( keys %clus_hits != grep { $_->[5] eq $clus_idx } @{$part} );

  my @clus_keys  = keys %clus_hits;
  my $clus_frags = list2frags( \@clus_keys );
  map { $clus_hits{ $_->{KEY} }->{FRAG} = $_ } @{$clus_frags};
  my @orig_clus = unique( map { $clus_hits{$_}->{PART}->[4] } keys %clus_hits );
  my $clus_dir = "RESULTS/$clus_idx";
  ## read in model ids of a final(merged) cluster, could be >1 orig clusters in case of merging
  my %model_ids = ();

  foreach my $f (sort(@modTreeFiles)) {

    my @model_fa = read_fasta_file("$f");
    map { $model_ids{$_} = 1 } @{ $model_fa[1] };
  }

  ## annotate %clus_hits with TYPE=MODEL or TYPE=BLASTCLUST
  my ( $model_map, $bc_map ) = getHitMap( $in_root_dir, \%clus_hits, \%model_ids );
  ## write out current cluster as partition
  open( PART, ">$clus_dir.cluster.part" );
  map { print PART join( " ", @{ $clus_hits{$_}->{PART} } ) . "\n"; } keys %clus_hits;
  close(PART);
  makeCol("$clus_dir.cluster.part");

  ############################################################################
  ## write cluster.all file with detailed cluster infos
  open( OUT, ">$clus_dir.cluster.all" );
  ## write model ids to cluster file
  foreach my $key ( sort { $clus_hits{$b}->{PART}->[1] <=> $clus_hits{$a}->{PART}->[1] } grep { $clus_hits{$_}->{TYPE} eq "MODEL" } keys %clus_hits ) {
    my $score = $clus_hits{$key}->{PART}->[1];
    my $clus  = $clus_hits{$key}->{PART}->[4];
    print OUT "CLUSTER  $clus_idx  CM_SCORE $score MODEL $clus ";
    my  $str = $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} };
    my @orId = split / /, $str;
    print OUT $orId[3] .  "\n";

  }

  ## write all other hit seqs
  foreach my $key ( sort { $clus_hits{$b}->{PART}->[1] <=> $clus_hits{$a}->{PART}->[1] } grep { $clus_hits{$_}->{TYPE} ne "MODEL" } keys %clus_hits ) {
    my $score = $clus_hits{$key}->{PART}->[1];
    my $clus  = $clus_hits{$key}->{PART}->[4];
    print OUT "CLUSTER $clus_idx   CM_SCORE $score " . $clus_hits{$key}->{TYPE} . " $clus ";
    my  $str = $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID}};
    my @orId = split / /, $str;
    print OUT $orId[3] .  "\n";
  }

  close(OUT);
  makeCol("$clus_dir.cluster.all");
  ############################################################################
  ## write fasta file for all hits in cluster
  my $fa_file_all = "$clus_dir.cluster.all.fa";
  writeClusterFasta( \@clus_keys, \%clus_hits, \@fa_scan, $fa_file_all );
  ############################################################################
  ## write stats file for global results stats
  my @orig_seqs = ();
  foreach my $key ( grep { $clus_hits{$_}->{PART}->[1] > 0 } sort { $clus_hits{$b}->{PART}->[1] <=> $clus_hits{$a}->{PART}->[1] } keys %clus_hits ) {
    my $header = $fa_scan[2]->{ $clus_hits{$key}->{FRAG}->{SEQID} };
    $header =~ /ORIGID\s+(\S+)/;
    push( @orig_seqs, $1 );
  }
  my @ids_unique = unique(@orig_seqs);
  open( STAT, ">$clus_dir.cluster.stats" );
  print STAT "CLUSTER $clus_idx SEQS " . scalar(@clus_keys) . " ";
  print STAT "IDS_UNIQUE " . scalar(@ids_unique) . " MODELS " . scalar(@orig_clus) . " ";
  # print STAT "IDS_UNIQUE_LIST ".join(";",@ids_unique); ## too long for column
  print STAT "\n";
  close(STAT);

}    ## for all @res_todo

exit;

#################################subs#################################
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
  ## check overlaps in more efficient way than all-vs-all
  ## achive this by special sorting of frags
  foreach my $f1 ( 0 .. $#{$fh1} ) {

    ## check if current (implies all) fragments in fh2 are "greater" that current f1
    next if ( $f2_curr_start_idx >= @{$fh2} || $fh1->[$f1]->{SEQID} lt $fh2->[$f2_curr_start_idx]->{SEQID} );

    ## get all fragments on the same seq in set $fh2
    if ( $fh1->[$f1]->{SEQID} ge $fh2->[$f2_curr_start_idx]->{SEQID} ) {
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
    }
    ## check all pairs, sorting allows to stop if f2 frag starts after frag f1
    foreach my $f2 (@f2_keys) {
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

sub makeCol {
  my $col_file = $_[0];
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

sub getTrueLocation {
  my $loc = $_[0];    ## data.locations loc identified by graphFasta.pl
  my $seq = $_[1];    ## seq of data.scan.fasta to get available length
  my $hit = $_[2];    ## fragment datastructure of cm hit

  return $loc if ( $loc =~ /MISS/ );

  ##hg19.chr1:949858-949920:+
  my @ent = split( ":", $loc );
  my @pos = split( "-", $ent[1] );
  ## GraphClust locations always have smaller start than end, but check anyway
  my $loc_len = abs( $pos[1] - $pos[0] );
  return $loc if ( $loc_len != length($seq) );
  my $new_strand = "";
  ## + + or - -
  $new_strand = "+" if ( $hit->{STRAND} eq $ent[2] );
  ## + - or - +
  $new_strand = "-" if ( $hit->{STRAND} ne $ent[2] );

  ## example for all 4 cases:
  ##
  ## chr1:100-250:-
  ## hit : 33- 87:+
  ## true:163-217:-

  ## chr1:100-250:+
  ## hit : 33- 87:-
  ## true:133-187:-

  ## chr1:100-250:+
  ## hit : 33- 87:+
  ## true:133-187:+

  ## chr1:100-250:-
  ## hit : 33- 87:-
  ## true:163-217:+

  ## overwrite old start and end
  if ( $ent[2] eq "+" ) {
    $pos[0] = $pos[0] + $hit->{START} - 1;
    $pos[1] = $pos[0] + $hit->{STOP} - $hit->{START};
  } else {
    $pos[0] = $pos[1] - $hit->{END} - 1;
    $pos[1] = $pos[0] + $hit->{STOP} - $hit->{START};
  }

  $loc = $ent[0] . ":" . $pos[0] . "-" . $pos[1] . ":" . $new_strand;

  return $loc;
}

sub writeClusterFasta {
  my $hits         = $_[0];
  my $clusHits     = $_[1];
  my $faScan       = $_[2];
  my $out_filename = $_[3];

  open( FA, ">$out_filename" );

  my %locations = ();
  open( LOCS, "FASTA/data.locations" );
  while ( my $line = <LOCS> ) {
    chomp($line);
    my @ent = split( " ", $line );
    $locations{ $ent[0] } = $ent[1];
  }
  close(LOCS);

  foreach my $frag ( @{$hits} ) {

    my $seq;
    my $start = $clusHits->{$frag}->{FRAG}->{START} - 1;
    my $len = $clusHits->{$frag}->{FRAG}->{STOP} - $clusHits->{$frag}->{FRAG}->{START} + 1;

    my $seqid = $clusHits->{$frag}->{FRAG}->{SEQID};
    if ( $clusHits->{$frag}->{FRAG}->{STRAND} eq "+" ) {
      $seq = substr( $faScan->[0]->{$seqid}, $start, $len );
    } elsif ( $clusHits->{$frag}->{FRAG}->{STRAND} eq "-" ) {
      $seq = substr( $faScan->[0]->{$seqid}, $start, $len );
      $seq =~ tr/AUGC/UACG/;
      $seq = reverse($seq);
    } else {
      die "Strand Error for frag $frag! Exit...\n\n";
    }

    my $part = $clusHits->{$frag}->{PART};
    my $fr   = $frag;
    my $loc  = $locations{ $clusHits->{$frag}->{FRAG}->{SEQID} };

    $loc = getTrueLocation( $loc, $faScan->[0]->{$seqid}, $clusHits->{$frag}->{FRAG} );

    $fr =~ s/#/_/g;
    print FA ">$fr" . " RESULT " . $part->[5] . " SCORE " . $part->[1] . " EVALUE " . $part->[2];
    print FA " CLUSTER " . $part->[4] . " LOC $loc ";

    ## add evaluation info if avaliable
    if ( @{$part} >= 12 ) {
      print FA "CLUS_CLASS " . $part->[7] . " CLASS " . $part->[6] . " CLASS_NAME " . $part->[11] . " CLASS_KEY " . $part->[8] . " CLASS_OL " . $part->[9] . " ";
    }

    print FA "$faScan->[2]->{$seqid}\n";
    print FA "$seq\n";
  }
  close(FA);
}

sub getHitMap {
  my $rootDir    = $_[0];
  my $clus_href  = $_[1];
  my $model_href = $_[2];
  opendir( BLASTCL, "FASTA/" );
  my @blast_cl = readdir(BLASTCL);
  closedir(BLASTCL);

  my %blast_cl = ();

  foreach my $file (@blast_cl) {
    next if ( $file !~ /intermediate.blast_clusters.(\d+)/ );
    open( IN, "$rootDir/FASTA/$file" );
    while ( my $line = <IN> ) {
      chomp($line);
      my @ent = split( " ", $line );

      foreach my $idx ( 0 .. $#ent ) {
        my @push = grep { $_ =~ /SEQ\d+#\d+#\d+#/ } map { $ent[$_] if ( $_ != $idx ) } 0 .. $#ent;

        #print "push($ent[$idx])".join(":",@push)."\n";
        if ( exists( $blast_cl{ $ent[$idx] } ) ) {
          push( @{ $blast_cl{ $ent[$idx] } }, @push );
        } else {
          $blast_cl{ $ent[$idx] } = \@push;
        }
      }

    }

    close(IN);
  }

  my @bc_keys  = keys %blast_cl;
  my $bc_frags = list2frags( \@bc_keys );

  my $hit_frags = [];
  map { push( @{$hit_frags}, $clus_href->{$_}->{FRAG} ) } keys %{$clus_href};
  @{$hit_frags} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$hit_frags};
  my $bc_ols = fragment_overlap( $hit_frags, $bc_frags, 0.51 );
  my %blast_hits = ();
  foreach my $ol ( @{$bc_ols} ) {
    $blast_hits{ $hit_frags->[ $ol->[0] ]->{KEY} } = 1;
    $clus_href->{ $hit_frags->[ $ol->[0] ]->{KEY} }->{TYPE} = "BLASTCLUST";
  }

  ###################################
  my $all_frags = read_fragments("FASTA/data.map");
  my %model_map   = ();
  my @model_frags = ();

  map { $model_map{ $_->{VALUE} }->{FRAG} = $_ if ( exists $model_href->{ $_->{VALUE} } ) } @{$all_frags};
  map { push( @model_frags, $model_map{$_}->{FRAG} ) } keys %model_map;

  @model_frags = sort { $a->{SEQID} cmp $b->{SEQID} } @model_frags;

  my $model_ols = fragment_overlap( \@model_frags, $hit_frags, 0.51 );

  foreach my $ol ( @{$model_ols} ) {
    $clus_href->{ $hit_frags->[ $ol->[1] ]->{KEY} }->{TYPE} = "MODEL";
    $model_map{ $model_frags[ $ol->[0] ]->{VALUE} }->{HIT} = $hit_frags->[ $ol->[1] ]->{KEY};
  }

  return ( \%model_map, \%blast_hits );
}

sub fasta2BED {

  my $fa_file  = $_[0];
  my $bed_file = $_[1];
  my $info     = $_[2];
  my @fa = read_fasta_file($fa_file);
  open( BED, ">$bed_file" );
  foreach my $id ( @{ $fa[1] } ) {
    my $header = $fa[2]->{$id};

    print $header. "\n";
    $header =~ /LOC\s+(\S+)/;
    my $loc = $1;

    if ( $loc =~ /MISS/ ) {
      next;
    }

    print "loc:$loc\n";
    my @ent = split( ":",  $loc );
    my @chr = split( /\./, $ent[0] );
    print "ent:" . join( "#", @ent ) . "#\n";
    print "chr:" . join( "#", @chr ) . "#\n";
    my @pos = split( "-", $ent[1] );
    my $strand = $ent[2];

    print BED $chr[1] . "\t" . $pos[0] . "\t" . $pos[1] . "\t$info\_$id\t0\t$strand\n";

  }
  close(BED);

  makeCol($bed_file);

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
