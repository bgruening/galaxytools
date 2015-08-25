#!/usr/bin/perl -w
use strict;




my $largest = 0;
my $contig = '';


if (@ARGV != 1) {
    print "$0 454Scaffolds.txt\n" ;
	exit ;
}

my $filenameA = $ARGV[0];

open (IN, "$filenameA") or die "oops!\n" ;
open OUT_1, ">", "$filenameA.scaffolds" or die "oops\n" ;
open OUT_2, ">", "$filenameA.singlet" or die "oops\n" ;


my %scaffold_lines = () ;
my %contigs_in_scaffold = () ;


while (<IN>) {


    my @r = split /\s+/, $_ ;
    
    next if $r[4] eq 'N' ;

    $scaffold_lines{$r[0]} .= "$r[5]\t$r[0]\n" ;
    $contigs_in_scaffold{$r[0]} = $r[3] ;
    


}


my $singlet = 0 ;


for my $scaffold (sort keys %contigs_in_scaffold ) {

    if ( $contigs_in_scaffold{$scaffold} == 1 ) {
	print OUT_2 "$scaffold_lines{$scaffold}" ;
	$singlet++ ;
    }
    else {
	print OUT_1 "$scaffold_lines{$scaffold}" ;
    }


}


print "a total of $singlet singlet scaffolds\n" ;
