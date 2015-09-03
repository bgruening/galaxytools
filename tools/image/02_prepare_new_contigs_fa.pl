#!/usr/bin/perl -w
#
#	30_walk_stage1.pl
#	Local assemblies of velvet
#
#
#
#


use strict;


if (@ARGV != 2) {
	print "02_prepare_new_contigs_fa.pl read_placed contigs.fa\n (will produce contigs.fa.new contigs.fa.new.list)\n\n" ;	
	exit ;
}


# sorting out different versions
my $read_placed = shift @ARGV ;
my $contig_fa = shift @ARGV ;

#my $read_placed = shift @ARGV ;


####################################################################################


# sort the contig and supercontig

open( IN, "$read_placed" ) or die "Cannot open $read_placed\n";

my %contig_in_supercontigs = () ;
my %supercontigs = () ;
my @supercontigs_ordered = () ;
my @contigs_ordered = () ;
my %contigs_names = () ;

while (<IN>) {
    chomp ;
    tr/\|/\./ ;


	next if /xxxxxx/ ;

	if (/(\S+)\s+(\S+)/) {	
		#print "$1\t$2\n" ;
		push @{ $contig_in_supercontigs{$2} } , "$1" ;
		push(@contigs_ordered, $1) ;

		unless ( $supercontigs{$2} ) {
			push(@supercontigs_ordered, $2) ;

		}	
		$supercontigs{$2}++ ;


		my $contig = $1 ;
		if ($contig =~ /(\d{5})\.(\S+)/ ) {
			
			$contigs_names{$2} = $contig ;
		}
		


	}
}
close(IN) ;



print "1. reading original contig.fa file!...\n" ;
open (IN, "$contig_fa") or die "oops! can not open $contig_fa\n" ;
my $seq_name = $1 ;
my $seq = '' ;
my %contig_seq = () ;

while (<IN>) {
    tr/\|/\./ ;



    if (/^>(\S+)/) {
	$seq_name = $1 ;

	while (<IN>) {
	    chomp ;
	    tr/\|/\./ ;



	    if (/^>(\S+)/) {
		#print "$seq_name\n" ;

		$contig_seq{$seq_name} = $seq ;
		$seq_name = $1 ;
		$seq = '' ;
	    }
	    else {
		$seq .= $_ ;
	    }

	}

    }
}

$contig_seq{$seq_name} = $seq ;
close(IN) ;



open OUT, ">", "image.contigs.fa" or die "Cannot open contigs.fa.new\n" ;
open OUT2, ">", "image.contigs.fa.list" or die "Cannot open contigs.fa.new.list\n" ;


for (keys %contig_seq ) {
	
	print OUT2 "$_\t$contigs_names{$_}\n" ;

	print OUT ">$contigs_names{$_}\n" ;
	print OUT "$contig_seq{$_}\n" ;


}


print "done!\ncontigs.fa.new and contigs.fa.new.list produced\n" ;
