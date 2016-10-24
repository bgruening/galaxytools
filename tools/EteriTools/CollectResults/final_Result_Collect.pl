#!/usr/bin/perl -w
use warnings;
use strict;
use Getopt::Long;
use Cwd qw(abs_path getcwd);
use Array::Utils qw(:all);
use POSIX qw(ceil floor);
use List::Util qw/ min max /;

print "stage 9: Final RESULTS - collect all found clusters";

my $RESULTS_DIR = "RESULTS";
system("mkdir $RESULTS_DIR");

my ($part_file) = @ARGV;

#my $part_file = "$RESULTS_DIR/partitions/final_partition.soft";

  if (!-e $part_file){
    print "No partition file found! Maybe 0 clusters found so far!";
    return;
  }

  ## part: array of [528,42.40,1.1,1.2,1,63,63]
  my $part = read_partition($part_file);
  my @res_clus_sort = sort { $a <=> $b } unique( map { $_->[5] } @{$part} );

  return if ( !@res_clus_sort );    ## do we have results at all

#  if ( !-e "$RESULTS_DIR/results.DONE" ) {

    system("perl gc_results_cluster.pl $part_file");


#  }##end if ( !-e "$RESULTS_DIR/results.DONE" )

  my $stats_file = "$RESULTS_DIR/cluster.final.stats";
  system("rm -f $stats_file");
  foreach my $clus (@res_clus_sort) {
    system("cat $RESULTS_DIR/$clus/cluster.stats >> $stats_file");
  }

  makeCol($stats_file);

  print "Results done!";

  #print "\nAll clusters can be found in:\n $ROOT_NAME/CLUSTERS\n\n";
#  print "To re-compute results please delete first file:\n $ROOT_NAME/RESULTS/results.DONE\n\n";
#  my $CMD_res = "MASTER_GraphCluster.pl --root $ROOT_NAME --results ";
#  $CMD_res .= "--sge " if ($in_USE_SGE);
#  print "Then invoke pipeline again with:\n $CMD_res\n";



##subs

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
