#!/usr/bin/env perl
use strict;

my $file = shift;
my $out1 = shift;
my $out2 = shift;
open FILE, "< $file" or die "can t open $file\n";
open OUT1, "> $out1" or die "can t write to $out1\n";
open OUT2, "> $out2" or die "can t write to $out2\n";

my $clusterNo = 1;
my @entry = ();
my %hash = ();

while(<FILE>){
        chomp $_;
	push @entry, $_;
	
	if (scalar(@entry) == 2){
		my ($id, $seq) = @entry;
		@entry = ();
		
		$hash{$seq} .= $id."\t";
	}	
}

foreach my $KEY (keys %hash){
	print OUT1 ">cluster$clusterNo\n$KEY\n";
	my @LINE = split/\t/,$hash{$KEY}; 
	foreach (@LINE){
		my @ID = split/>/,$_;
		print OUT2 ">cluster$clusterNo".":$ID[1]\n$KEY\n";
	}
	$clusterNo++;	
}

