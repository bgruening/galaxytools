#!/usr/bin/env perl

use strict;

##------------------------------------------------------------
## options
use Getopt::Long;
use Pod::Usage;

my $help;
my $man;
my $quiet;
my $verbose;

my $quantile = 0.9993;    #original RNAclust

#my $quantile=1;

## Getopt::Long::Configure("no_ignore_case");

GetOptions(
    "verbose"    => \$verbose,
    "quiet"      => \$quiet,
    "help"       => \$help,
    "man"        => \$man,
    "quantile=f" => \$quantile
) || pod2usage(2);

pod2usage(1) if $help;
pod2usage( -exitstatus => 0, -verbose => 2 ) if $man;

## ------------------------------------------------------------
## main part

my @lines  = <STDIN>;
my @scores = ();

for my $line (@lines) {

    #if($line =~/(\d+)$/) {
    #dont push lines which have inf, as score (necessary for local clustering
    if ( ( $line =~ /^\S+\s+\S+\s+(\S+)\s+/ ) && ( $line !~ /inf/ ) ) {
        push @scores, $1;
    }
}

my @sscores = sort { $a <=> $b } @scores;

my $m = $sscores[ $quantile * ($#scores) ];

for my $line (@lines) {

    #if($line =~/(.*)\s(\d+)$/) {
    if ( $line =~ /^(\S+\s+\S+)\s+(\S+)/ && ( $line !~ /inf/ ) ) {
        print "$1 " . ( ( $m - $2 > 0 ) ? ( $m - $2 ) : 0 ) . "\n";
    }
    else {
        print $line;
    }
}

## ------------------------------------------------------------

__END__

=head1 NAME

rnaclustScores2Dist.pl

=head1 SYNOPSIS

rnaclustScores2Dist.pl [options]

Options:

        --help           brief help message

        --man            full documentation

        --verbose        be verbose

        --quiet          be quiet

        --quantile       the quantile

=head1 DESCRIPTION

Read list of scores from stdin and writes list of distances to stdout

=cut
