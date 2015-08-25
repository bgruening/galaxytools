#!/usr/bin/perl
#
#	to incorporate the information from velvet contig (in .afg)
#		to original ace file (which has been converted in .afg)
#
#	author: jit
#
#
#
#


use strict;



if (@ARGV != 2) {
	print "xxx.pl walk2 folder fastas outprefix\n\n" ;
	exit ;
}



my $folder = shift @ARGV ;
my $prefix = shift @ARGV ;


open my $IN, q[-|] , qq[ cat $folder/walk2* ];
open my $fasta , ">" , "$prefix.fa" or die "oops\n" ;
open my $out , ">" , "$prefix.out" or die "oops\n" ;


while (<$IN>) {
	chomp ;

	my $vel_contig ;
	if (/(NODE_\S+)\s+\d+\s+2/) {
		$vel_contig = $1 ;
	}
	else {
		next ;
	}

	print $out "Found gap closed!\n" ;

	my $resultline = <$IN> ;
	my @result = split /cigar::50/ , $resultline ;
	shift(@result) ;

	my $count = 0 ;
	my $fasta_left ;
	my $fasta_right ;

	foreach (@result) {
		print $out "$_\n" ;
		my @tmp = split /\s+/ ;

		if ($count == 0 ) {
			$fasta_left = $tmp[5] ;
		}
		else {
			$fasta_right = $tmp[5] ;
		}

		$count++ ;
	}



	print $fasta ">$fasta_left.$fasta_right.$vel_contig\n" ;

	#parse out the fasta
	#command is like this...
	#fasta2singleLine.pl ite31_again1/velvet_31.25.10/contig_1213.contig_1214/newcontigs.fa | grep '$vel_contig' -A 1 | grep '$vel_contig' -v

	print $fasta `fasta2singleLine.pl $folder/velvet_*/$fasta_left.$fasta_right/newcontigs.fa | grep '$vel_contig' -A 1 | grep '$vel_contig' -v` ;
		


}



print "all done!\n" ;

