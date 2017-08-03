#!/usr/bin/perl

use strict;
use warnings;
#se Data::Dumper;

my $in = shift;
my $out = shift;

open IN, "< $in" or die "can t open $in\n";
open OUT, "> $out" or die "can t open $out\n";


my @ids  = ();
my %seqs = ();
my $id   = '';

while(<IN>){
	chomp;
  	if ($_ =~ m/^>/){
    		$id = $_;
    		push @ids, $id;
 	}
 	else{
    		$seqs{$id}.=$_;
	}
}

foreach my $id (@ids){
  	if($id !~ m/pseudo/i){
		print OUT "$id\n";
  		print OUT join("", lc($seqs{$id}), "cca")."\n";
 	}
}
