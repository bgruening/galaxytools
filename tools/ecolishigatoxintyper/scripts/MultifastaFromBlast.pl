#!/usr/bin/env perl
use strict;
use warnings;
use English;

# Parse arguments
my ($inputs, $output) = @ARGV;
# Run program
unlink $output;
my @infiles = split( /,/, $inputs );
foreach(@infiles) {
  transformFileContent($_, $output);
}

# read input file and write multifasta with sequence (forward or reverse complement)
sub transformFileContent {
  my ($infile, $outfile) = @_;
  open my $if, '<', $infile or die "Cannot open : $infile!";
  open my $of, '>>', $outfile or die "Cannot open : $outfile!";
  my @lines = <$if>;
  close $if;
  foreach(@lines) {
    my @elems = split( /\t/, $_ );
    print $of ">$elems[0]\n";
    chomp $elems[2];
    if ($elems[1] == 1) {
      print $of "$elems[2]\n";
    }
    else {
	  my $revcomp = reverseComplement($elems[2]);
      print $of "$revcomp\n";
    }
  }
  close $of;
  return 0;
}

# calculate reverse complement
sub reverseComplement {
  my ($DNA) = @_;
  my $revcom = reverse $DNA;
  $revcom =~ tr/ACGTacgt/TGCAtgca/;
  return $revcom;
}
