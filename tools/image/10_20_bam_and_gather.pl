#! /usr/bin/perl -w
#
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: JIT
# Description: 
#		a script that map the Solexa reads back to the reference/contigs
#		and parition them based on either ends of the contigs
#
#
#



use strict;
use warnings;



if (@ARGV != 4 ) {
	print "10_bam_and_gather.pl seqcoords bam insert_size unmapped_only[1]|or_not[0]\n\n" ;
	exit;
}


my $seqcoord = shift ;
my $output = shift ;
my $range = shift ;
my $unmapped = shift ;


my %BINS = () ;

#chr19 sequence coordinates
#mouse_ncbim37_19.0 3000001 6685768 , len = 3685767
#mouse_ncbim37_19.1 6685869 6688105 , len = 2237
#mouse_ncbim37_19.2 6888106 6904236 , len = 16130
#mouse_ncbim37_19.3 6904337 61342430, len = 54438093


open( IN, "$seqcoord" ) or die "Cannot open output: $seqcoord \n";
while (<IN>) {
	chomp;
	
	my @line = split /\s+/ ;
	$BINS{$line[1]}{"$line[0].L"} = "$line[2]." . ( $line[2] + $range ) ;
	$BINS{$line[1]}{"$line[0].R"} = ( $line[3] - $range ) . ".$line[3]"   ;

}
close(IN) ;





# start parsing..
if ( $unmapped == 0 ) {
	open( IN, "samtools view $output |" ) or die "Cannot open piped output of samtools view $output \n";
}
else {
	open( IN, "samtools view -F 0x2 $output |" ) or die "Cannot open piped output of samtools view $output \n";
}
#print "parsing output....\n" ;




my %contig_reads ;
my $total_reads = 0 ;

#my $count = 0 ;


while (<IN>) {

	if (/^\@SQ/ ) {
		next ;
	}

	chomp;
	my @sam = split /\s/ , $_ ;

	#print "@sam\n" ;


	# exclude some of the dodgy reads
	if ( $sam[1] == 113 || $sam[1] == 177 ) {
		next ;
	}
	elsif ( $sam[1] == 161 || $sam[1] == 81 ) {
		next ;
	}
	elsif ( $sam[1] == 97 || $sam[1] == 145 ) {
		next ;
	}
	elsif ( $sam[1] == 65 || $sam[1] == 129 ) {
		next ;
	}
	




	my $isBin = locate_bin($sam[2],$sam[3],$sam[7]) ;
	if ($isBin) {
		print "$_ $isBin\n" ;
	}


	#last ;


}



sub locate_bin {


my $chr = shift ;
my $coord_L = shift ;
my $coord_R = shift ;

	for my $bin (keys %{ $BINS{$chr} }) {
		my @pos = split /\./ , $BINS{$chr}{$bin}  ;
		
		#print "$coord $pos[0] $pos[1]\n" ;

		if ( $coord_L >= $pos[0] && $coord_L <= $pos[1] ) {
			return $bin ;
		}
		elsif ( $coord_R >= $pos[0] && $coord_R <= $pos[1] ) {
			return $bin ;
		}
	}


}




