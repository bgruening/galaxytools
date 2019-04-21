#!/usr/bin/perl -w

use strict;
use Getopt::Long;
use POSIX qw(ceil floor);

my $maxLength     = 400;
my $desiredLength = 200;
my $help          = 0;

GetOptions(
  "maxLength:i"     => \$maxLength,
  "desiredLength:i" => \$desiredLength,
  "help"            => \$help,
  "h"               => \$help
);

if ($help) {
  print STDERR "\nbreakMAF [options] < input.maf > output.maf\n\n";
  print STDERR "Options:\n";
  print STDERR "  maxLength:     Break all blocks longer than that (default: 400 columns)\n";
  print STDERR "  desiredLength: Try to create blocks of this size (default: 200 columns)\n";
  exit(1);
}

$/ = 'a score=';

while ( my $block = <> ) {

  $block = 'a score=' . $block;

  my $currAln = readMAF($block)->[0];

  next if not defined($currAln);

  my $length = length( $currAln->[0]->{seq} );

  if ( $length > $maxLength ) {

    my $breakN      = int( $length / $desiredLength );
    my $chunkLength = ceil( $length / $breakN );

    my $from = 0;

    while (1) {
      my $to = $from + $chunkLength;

      if ( $to > $length ) {
        $to = $length;
      }

      my $slice = sliceAlnByColumn( $currAln, $from, $to );

      print formatAln( $slice, 'maf' ), "\n";

      $from = $to;

      last if $from == $length;
    }

  } else {
    print formatAln( $currAln, 'maf' ), "\n";
  }
}

######################################################################
#
# readMAF($string)
#
# Converts the MAF in string to internal alignment format. Returns
# list of usual array references.
#
######################################################################

sub readMAF {

  my $string = shift;

  return [] if $string eq '';

  my @input = split( "\n", $string );

  #open(FILE,"<$file") || die("Could not read $file ($!)");

  my @outAlns = ();
  my @aln     = ();

  foreach my $i ( 0 .. $#input ) {

    $_ = $input[$i];

    next if (/^\s?\#/);
    next if (/^\s?a/);

    if (/^\s?s/) {
      ( my $dummy, my $name, my $start, my $length, my $strand, my $fullLength, my $seq ) = split;

      my $end = $start + $length;

      ( my $org, my $chrom ) = split( /\./, $name );

      push @aln, {
        name   => $name,
        org    => $org,
        chrom  => $chrom,
        start  => $start,
        end    => $end,
        seq    => $seq,
        strand => $strand
        };
    }

    if ( /^\s?$/ and @aln ) {
      push @outAlns, [@aln];
      @aln = ();
    }
    if ( ( not defined $input[ $i + 1 ] ) and (@aln) ) {
      push @outAlns, [@aln];
      @aln = ();
    }
  }
  return \@outAlns;
}

sub formatAln {

  my @aln = @{ $_[0] };

  my $format = $_[1];

  $format = lc($format);

  my @alnSeqs  = ();
  my @alnNames = ();

  my $counter = 1;

  foreach my $row (@aln) {

    my $name = "seq$counter";
    $counter++;
    if ( $row->{name} ) {
      $name = ( $row->{name} );
    } elsif ( $row->{org} and $row->{chrom} ) {
      $name = "$row->{org}.$row->{chrom}";
    }

    my $start  = $row->{start};
    my $end    = $row->{end};
    my $strand = $row->{strand};

    my $pos = '';

    if (  defined $start
      and defined $end
      and defined $strand
      and ( $format ne 'phylip' ) ) {
      ( $start, $end ) = ( $end, $start ) if ( $strand eq '-' );
      $pos = "/$start-$end";
    }

    push @alnNames, "$name$pos";
    push @alnSeqs,  $row->{seq};

  }

  my $output = '';

  if ( $format eq 'clustal' ) {

    $output = "CLUSTAL W(1.81) multiple sequence alignment\n\n\n";
    my $maxName = 0;

    foreach my $name (@alnNames) {
      $maxName = ( $maxName < length($name) ) ? length($name) : $maxName;
    }

    for my $i ( 0 .. $#alnNames ) {
      my $buffer = " " x ( ( $maxName + 6 ) - length( $alnNames[$i] ) );
      $alnNames[$i] .= $buffer;
    }
    my $columnWidth = 60;
    my $currPos     = 0;
    my $length      = length( $alnSeqs[0] );

    while ( $currPos < $length ) {
      for my $i ( 0 .. $#alnNames ) {
        $output .= $alnNames[$i];
        $output .= substr( $alnSeqs[$i], $currPos, $columnWidth );
        $output .= "\n";
      }
      $output .= "\n\n";
      $currPos += $columnWidth;
    }
  } elsif ( $format eq 'fasta' ) {
    foreach my $i ( 0 .. $#alnNames ) {
      my $name = $alnNames[$i];
      my $seq  = $alnSeqs[$i];
      $seq =~ s/(.{60})/$1\n/g;
      $output .= ">$name\n$seq\n";
    }
  } elsif ( $format eq 'phylip' ) {
    my $length = length( $alnSeqs[0] );
    my $N      = @alnSeqs;
    $output .= "  $N $length\n";
    foreach my $i ( 0 .. $#alnNames ) {
      my $name = $alnNames[$i];
      my $seq  = $alnSeqs[$i];
      $seq =~ s/(.{60})/$1\n/g;
      $output .= "$name\n$seq\n";
    }
  } elsif ( $format eq 'maf' ) {
    $output .= "a score=0\n";
    foreach my $row (@aln) {
      my $length = $row->{end} - $row->{start};
      $output .= "s $row->{org}.$row->{chrom} $row->{start} $length $row->{strand} 0 $row->{seq}\n";
    }
  }
  return $output;

}

######################################################################
#
# sliceAlnByColumn(\@aln ref-to-alignment, $start int, $end int)
#
# Returns slice of alignment specified by alignment column.
#
# \@aln ... alignment in list of hash format
# $start, $end ... slice to cut
#
# Returns reference to alignment in list of hash format. This is a new
# alignment, i.e. the input is not sliced in place
#
######################################################################

sub sliceAlnByColumn {

  my @aln = @{ $_[0] };
  shift;
  ( my $start, my $end ) = @_;

  # correct ends without warning if outside of valid range
  $start = 0 if ( $start < 0 );
  $end = length( $aln[0]->{seq} ) if ( $end > length( $aln[0]->{seq} ) );

  #my @newAln=@aln;

  # make deep copy of list of hash
  my @newAln = ();
  foreach (@aln) {
    push @newAln, { %{$_} };
  }

  foreach my $i ( 0 .. $#newAln ) {

    my $oldStart = $newAln[$i]->{start};
    my $oldEnd   = $newAln[$i]->{end};

    $newAln[$i]->{start} = alnCol2genomePos( $newAln[$i]->{seq}, $oldStart, $start );
    $newAln[$i]->{end}   = alnCol2genomePos( $newAln[$i]->{seq}, $oldStart, $end - 1 ) + 1;
    $newAln[$i]->{seq} = substr( $newAln[$i]->{seq}, $start, $end - $start );

  }

  return ( [@newAln] );

}

######################################################################
#
# alnCol2genomePos($seq string, $start int, $col int)
#
# Calculates the genomic position corresponding to a column in an
# alignment.
#
# $seq ... sequence from alignment (i.e. letters with gaps)
# $start ... Genomic position of first letter in $seq
# $col ... column in the alignment that is to be mapped
#
# Returns genomic position. No error handling, so $col must be a valid
# column of the string $seq.
#
#######################################################################

sub alnCol2genomePos {

  ( my $seq, my $start, my $col ) = @_;

  $seq =~ s/\./-/g;    #Convert all gaps to "-"

  my $newPos = $start;

  # if gap only...
  return $start if ( $seq =~ /^-+$/ );

  ( my $tmp ) = $seq =~ /(-*)[^-]/;

  my $leadingGaps = length($tmp);

  # if column is in gap before first letter,
  # return position of the first letter
  return $start if ( $col < $leadingGaps );

  $newPos = $start - 1;

  for my $i ( $leadingGaps .. $col ) {
    $newPos++ if ( ( my $t = substr( $seq, $i, 1 ) ) ne '-' );
  }
  return $newPos;
}

