#!/usr/bin/perl -w
use strict;




my $largest = 0;
my $contig = '';


if (@ARGV != 1) {
	print "contig_length_fastq fasta\n\n" ;
	exit ;
}

my $filenameA = $ARGV[0];

open (IN, "$filenameA") or die "oops!\n" ;

	my $read_name = '' ;
	my $read_seq = '' ;

	while (<IN>) {
	    if (/^>(\S+)/) {
		$read_name = $1 ;
		$read_seq = "" ;
		
		while (<IN>) {

			if (/^>(\S+)/) {
			    
			    print "$read_name\t" . length($read_seq) . "\n" ;

			    $read_name = $1 ;
			    $read_seq = "" ;



			}
			else {
			    chomp ;
			    $read_seq .= $_ ;
			}


		}

	    }
	}

close(IN) ;

print "$read_name\t" . length($read_seq) . "\n" ;
