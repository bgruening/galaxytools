#!/usr/bin/perl -w
use strict;




my $largest = 0;
my $contig = '';


if (@ARGV < 4) {
	print "\n\nIMAGE: contigs to scaffolds script\n" ;
	print "Author: jit\@sanger.ac.uk\n\n" ;

	print "Usage: $0 new.fa new.read.placed gap_length min_length_contig_to_exclude output_prefix\n\n" ;
	print "\tDescription:\tconverts imaged contigs into scaffolds (joined by Ns)\n" ; 
	print "\t\t\tnew.fa and new.read.placed are files from the last iteration\n" ; 
	
	print "\n\tExample: $0 Ite_10/new.fa Ite_10/new.read.placed 300 500 scaffolds\n\n" ; 
	print "\t\tThis will create scaffolds from the Ite_10 folder where contigs less then 500bp will be excluded. \n\t\tThe contigs will be joined by gaps of 300 Ns\n" ; 
	print "\t\tA final file called: scaffolds.fa will be produced\n\n\n" ; 

	exit ;
}

my $filenameA = $ARGV[0];
my $filenameB = $ARGV[1];
my $gap_len =  $ARGV[2];
my $min_len = $ARGV[3] ;
my $output =  $ARGV[4];







my %contig_len = () ;
my %contig_seq = () ;
my %contig_present = () ;


#opening file handles
open (IN, "$filenameA") or die "oops!, can not open contig.fasta: $filenameA\n" ;


open my $out_fa, ">" , "$output.fa" or die "can not open $output.fa\n" ;
open my $out_pos, ">" , "$output.pos" or die "can not open $output.pos\n" ;


open my $old_read_placed, "< $filenameB" or die "can not open $filenameB!\n" ;
print "readplaced file is : $filenameB\n" ;
print "gaps will be removed in the consensus!!\n" ;

while (<IN>) {

    # print "$_" ;

    if (/^>(\S+)/) {

	my $seq_name = $1 ;
	my $seq = <IN> ;
	chomp($seq) ;
	my $seq_len = length($seq) ;       

	$seq_name =~ s/^\d\d\d\d\d\.// ;

	$seq =~ s/-//g ;

	$contig_len{$seq_name} = $seq_len ;
	$contig_seq{$seq_name} = $seq ;
	$contig_present{$seq_name}++ ;

    }
}
close(IN) ;

#print "\#\#the largest length is: $contig with $largest bp \n" ;


my %contigs_in_supercontigs ;
my %contigs_found ;
my %supercontigs_found ;
my @supercontigs_ordered ;



while(<$old_read_placed>) {
	next if /xxxxx/ ;

	s/^\d\d\d\d\d\.// ;
	
	my @read_placed_split = split /\s+/ , $_ ;
	
	my $contig_original_id = '' ;
	
	if ( $read_placed_split[0] =~ /\.contig_(\d+)/ ) {
	    $contig_original_id = $1 ;
	}
	
	
	
	if ( $contigs_found{$read_placed_split[0]} ) {
	    
	}
	else {
	    $contigs_in_supercontigs{$read_placed_split[1]} .= "$read_placed_split[0] " ;
	    $contigs_found{$read_placed_split[0]}++ ;
	}
	
	unless ( $supercontigs_found{$read_placed_split[1]} ) {
	    push(@supercontigs_ordered, $read_placed_split[1]) ;
	    $supercontigs_found{$read_placed_split[1]}++ ;
	}
	
	
    }
close($old_read_placed) ;


my $count = 0  ;


foreach my $supercontig (@supercontigs_ordered) {

		my @contigs = split /\s+/ ,  $contigs_in_supercontigs{$supercontig};

		#print "doing SC: $supercontig\n" ;

		#my $SC_number = 0 ;
		#if ($supercontig =~ /(\d+)/) {
		#	$SC_number = $1 ;
		#}

		#$SC_number = sprintf("%04s", $SC_number) ;
		my $seq = '' ;
		my $seq_present = 0 ;

		my $pos = 1 ;
		my $previous_contig = '' ;

		my $start = 0 ;

		for (my $i = 0 ; $i < @contigs ; $i++) {
			my $contig = $contigs[$i] ;
			$count++ ;
			
			
			
			if ( $contig_seq{$contig} && length($contig_seq{$contig}) >= $min_len ) {

			    $seq_present = 1 ;
			    
			    if ( $start == 1 ) {
				
				my $gap = 'N' x $gap_len ;
			
				
				$seq .= "$gap$contig_seq{$contig}" ;

				printf $out_pos "$supercontig\tGap\t$pos\t" . ($pos+length($gap)-1) . "\n" ;								    
				$pos += length($gap) ;
				
				printf $out_pos "$supercontig\t$contig\t$pos\t" . ($pos+length($contig_seq{$contig})-1) . "\n" ;
				$pos += length($contig_seq{$contig}) ;
				
				
			    }
			    else {
				$seq .= $contig_seq{$contig} ;
				
				
				printf $out_pos "$supercontig\t$contig\t$pos\t" . ($pos+length($contig_seq{$contig})-1) . "\n" ;
				$pos += length($contig_seq{$contig}) ;
				
				$start = 1 ;
				
				
			    }
			    
			}
			else {

			    print "$contig leng: " . length($contig_seq{$contig}) . "\n ignored\n" ;


			}
			


		}

		if ( $seq_present == 1 ) {
		    $supercontig =~ s/\|/\./gi ; 

		    my $seq_len = length($seq) ; 
		    $supercontig =~ s/size(\d+)/size$seq_len/gi ; 



			print $out_fa ">$supercontig\n$seq\n" ;
		}




		#last;

}

print "a total of $count contigs have been used..\n" ;
print "all done!\n\n" ;






