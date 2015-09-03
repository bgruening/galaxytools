#!/usr/bin/perl -w
#
#	check Ns and change the name of the contig.fa
#	if the label says N, the length of the corresponding sequence is (N+k-1) 
#	the script will change the label to N+k-1
#
#


use strict;




if (@ARGV != 3) {
	print "check_velvet_contig.pl contigs.fa kmer newcontigs.fa\n\n" ;
	exit ;
}


my $contig_file = shift @ARGV ;
my $kmer = shift @ARGV ;
my $new_contig = shift @ARGV ;
my %contigs ;

open my $out ,">","$new_contig" or die "can't open output file!\n" ;


open( IN, "$contig_file" ) or die "Cannot open $contig_file\n";
my $seq_name = ''; my $seq = '';
my @seq_name_ordered ;
while (<IN>) {
	chomp;

	if (/>(\S+)/) {
		if ($seq_name ne '') {
			$contigs{$seq_name} = $seq ;
			$seq = '' ;
		}		
		$seq_name = $1 ;
		push(@seq_name_ordered, $seq_name) ;	
	}
	else {
		$seq .= $_ ;
	}
}
$contigs{$seq_name} = $seq ;


my $containN = 0 ;
foreach my $seq_names (@seq_name_ordered) {
	if ($contigs{$seq_names} =~ /N/ ) {
		print "$seq_names contains Ns!\n" ;
		$containN = 1 ;
	}
}


	foreach my $seq_names (@seq_name_ordered) {

		my $seq_names_changed = "" ;
		if ($seq_names =~ /(\S+)_length_(\d+)(\S+)/) {
				$seq_names_changed = "$1" . "_length_" . ($kmer+$2-1)  ;
		}

		if ($containN == 1 ) {
			my @sequences = split /N+/ , $contigs{$seq_names} ;
			my $count = 1 ;
			foreach (@sequences) {
				print $out ">$seq_names_changed.$count\n$_\n" ;
				$count++ ;
			}
		}
		else {
			print $out ">$seq_names_changed\n$contigs{$seq_names}\n" ;
		}	
	}
	close($out);



print "done done!\n" ;


