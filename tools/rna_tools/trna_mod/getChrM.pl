#!/usr/bin/perl
use strict;

my $file = shift;
my $out = shift;
open FILE, "< $file" or die "can t open $file\n";
open OUT, "> $out" or die "can t write to $out\n";


my $ID = '';
while(<FILE>){	
	if($_ =~m/^>/){
		$ID = $_;
	}
	if($ID =~m/chrM/i){
		print OUT $_;
	}	
}


