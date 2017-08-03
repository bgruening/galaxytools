#!/usr/bin/env perl

use strict;

my $file1 = shift; #pre-tRNAs
my $file2 = shift; #sam
my $out   = shift; #filtered sam
open FILE1, "< $file1" or die "can t open $file1\n";
open FILE2, "< $file2" or die "can t open $file2\n";
open OUT,   "> $out" or die "can t open $out\n";

my %ID   = ();
my %SEQ  = ();
my @SAM  = ();


while(<FILE1>){
        chomp $_;
	
	if($_=~ m/>/){
		my @id = split/>/,$_;
		$ID{$id[1]} = 0;
	}
}


while(<FILE2>){
        chomp $_;
	if($_=~ m/@/){
		print OUT $_."\n";
	}

	my @LINE = split/\t/,$_;

	if(!exists $SEQ{$LINE[0]}){
		$SEQ{$LINE[0]} = 0;
	}
	if(!exists $ID{$LINE[2]}){
		$SEQ{$LINE[0]} = 1;
	}
	push @SAM,$_;
}


foreach (@SAM){
	my @LINE = split/\t/,$_;
	if($SEQ{$LINE[0]} eq 0){
		print OUT $_."\n"; 
	}
}


