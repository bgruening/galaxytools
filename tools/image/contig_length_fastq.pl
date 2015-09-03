#!/usr/bin/perl -w
use strict;
use lib qw(/nfs/team81/jit/repository/pathogen/user/jit );
use strict;
use AssemblyTools;


my $largest = 0;
my $contig = '';


if (@ARGV != 1) {
	print "contig_length_fastq fastq\n\n" ;
	exit ;
}

my $filenameA = $ARGV[0];

open (IN, "$filenameA") or die "oops!\n" ;

while (<IN>) {

    # print "$_" ;

    if (/^\@(\S+)/) {

	my $seq_name = $1 ;
	my $seq = <IN> ;
	chomp($seq) ;
	my $seq_len = length($seq) ;       

	print "$seq_name\t$seq_len\n" ;

	if ($seq_len > $largest) {
	    $largest = $seq_len ;
	    $contig  = $seq_name;
	}

    }
	my $skip = <IN> ;
	$skip = <IN> ;

    
    #last;


}

print "\#\#the largest length is: $contig with $largest bp \n" ;
