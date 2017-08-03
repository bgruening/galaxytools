#!/usr/bin/perl
use strict;


my $in = shift;         
my $out = shift;        

open IN, "< $in" or die "can t open $in\n";
open OUT, "> $out" or die "can t open $out\n";


while(<IN>) {
        chomp $_;

	if($_ !~ m/pseudo/i){
		my @line = split(/\t/,$_);

		if($line[9] eq 2){
			my @SP = split(/,/,$line[10]);
			my $block = (($line[2] - $line[1]) - $SP[1]) +50;
			print OUT join("\t",$line[0],$line[1]-50,$line[2]+50,$line[3],$line[4],$line[5],$line[1]-50,$line[2]+50,$line[8],$line[9]);
			print OUT "\t";
			print OUT join (",",$SP[0]+50,$SP[1]+50);			
			print OUT "\t";
                	print OUT join (",",0,$block);
			print OUT "\n";
		}
		else{
			print OUT join("\t",$line[0],$line[1]-50,$line[2]+50,$line[3],$line[4],$line[5],$line[1]-50,$line[2]+50,$line[8],$line[9],$line[10]+100,$line[11]);
			print OUT "\n";
		}
	}
}
