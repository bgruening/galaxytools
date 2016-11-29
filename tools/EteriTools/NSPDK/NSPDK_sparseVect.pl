#!/usr/bin/perl
use List::Util qw/ min max /;
use POSIX qw(ceil floor);
#use Array::Utils qw(:all);

#my $CI = 1; ## iteration num
my $SVECTOR_DIR = "SVECTOR";


my ($data_fasta, $gspan, $rad, $dist, $path) = @ARGV;

#print $nspdk_knn_center;
my $group_gspan = "group.gspan";
#my $group_gspan = $gspan;
#print "gspan = $gspan \n";
#print system("file $gspan") ;
#print "\n";
#print system("file $data_fasta") ;

system("bzcat $gspan > $group_gspan");
#$group_gspan = $gspan;##################################

my @fa  = read_fasta_file($data_fasta);
my $num_seqs = @{ $fa[1] };
print "\nNumber of sequences in data.fasta: " . $num_seqs . "\n";

die "Fasta file $data_fasta contains only $num_seqs sequences! Exit...\n\n" if ( $num_seqs <= 2 );

system("$path./NSPDK -R $rad -D $dist -gt DIRECTED -fg $group_gspan -of") == 0 or die "nspdk chi ashxatum\n";
system("mkdir -p $SVECTOR_DIR");
system("cat $group_gspan.feature_bin > $SVECTOR_DIR/data.svector");
#system("cat $group_gspan.feature_bin > $SVECTOR_DIR/data.svector.$CI");
system("rm $group_gspan.feature_bin");
system("rm $group_gspan.feature");



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
