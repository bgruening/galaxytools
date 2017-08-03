#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;

my $in = shift;
my $out = shift;

open IN, "< $in" or die "can t open $in\n";
open OUT, "> $out" or die "can t open $out\n";



while(<IN>){
	chomp;
	if($_ =~ m/^>/){
 		my @LINE = split(/:/,$_);
  		print OUT $LINE[0]."\n";		
	}
	else{
		print OUT $_."\n";
	}
}

