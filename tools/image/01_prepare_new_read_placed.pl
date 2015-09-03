#!/usr/bin/perl -w
#
#	30_walk_stage1.pl
#	Local assemblies of velvet
#
#
#
#


use strict;


if (@ARGV != 1) {
	print "01_prepare_new_read_placed.pl read_placed \n (will produce read.placed.new)\n\n" ;	
	exit ;
}


# sorting out different versions
my $read_placed = shift @ARGV ;
#my $read_placed = shift @ARGV ;


####################################################################################


# sort the contig and supercontig

open( IN, "$read_placed" ) or die "Cannot open $read_placed\n";

my %contig_in_supercontigs = () ;
my %supercontigs = () ;
my @supercontigs_ordered = () ;

while (<IN>) {
    chomp ;
    tr/\|/\./ ; 


	if (/(\S+)\s+(\S+)/) {	
		#print "$1\t$2\n" ;
		push @{ $contig_in_supercontigs{$2} } , "$1" ;

		unless ( $supercontigs{$2} ) {
			push(@supercontigs_ordered, $2) ;
		}	

		$supercontigs{$2}++ ;

	}
}
close(IN) ;

open OUT, ">", "image.read.placed" or die "Cannot open read.placed.new\n" ;


#my $tmp_SC = 0 ;


foreach my $supercontig ( @supercontigs_ordered ) {

my @contigs = @{ $contig_in_supercontigs{$supercontig} } ;

print "doing SC: $supercontig\n" ;

my $count = 1 ;

	# for every contig pairs in supercontig
	for ( my $i = 0 ; $i < ( $#contigs + 1) ; $i++ ) {

		if ($i == 0 ) {
			print OUT "xxxxxx\t $supercontig\n" ;
		       
		}

		my $contig_count = sprintf("%05d", $count) ;
		print OUT "$contig_count.$contigs[$i]\t $supercontig\n" ;
		$count++ ;


		if ($i == $#contigs ) {	    
			print OUT "xxxxxx\t $supercontig\n" ;
		}

	
		
	
	}	


} #end of foreach loop





