#!/usr/bin/perl -w
use strict;




my $largest = 0;
my $contig = '';


my $path = '' ;
if ($0 =~ /(\S+)\/get_fasta_containNs_toIMAGE_format.pl/ ) {
    #print "The path is $1\n" ;                                                                                                                                                
    $path = $1 ;
}


if (@ARGV != 3) {
	print "$0 fasta threshold smallest_contig\n\n" ;
	exit ;
}

my $filenameA = $ARGV[0] ;
my $threshold = $ARGV[1] ;
my $smallest =  $ARGV[2] ;

print "checking all the paths..\n" ;
checkpath("$path/fasta2singleLine.pl") ;


#system("fasta2singleLine.pl $filenameA tmp.fa") ;
system("$path/fasta2singleLine.pl $filenameA tmp.fa") ;


open (IN, "tmp.fa") or die "oops!\n" ;
open OUT , ">", "$filenameA.above$threshold.fasta" or die "oops!\n" ;
open READ , ">", "$filenameA.above$threshold.read.placed" or die "oops!\n" ;

#open COORD , ">", "$filenameA.above$threshold.seqcoord" or die "oops!\n" ;
#open OUT2 , ">", "$filenameA.below$threshold.fasta" or die "oops!\n" ;
#open READ2 , ">", "$filenameA.below$threshold.read.placed" or die "oops!\n" ;



while (<IN>) {

    # print "$_" ;

    if (/^>(\S+)/) {

	my $seq_name = $1 ;
	my $seq = <IN> ;
	chomp($seq) ;
	my $seq_len = length($seq) ;       

	$seq = uc($seq) ; 

	if ($seq =~ /N/ ) {

		my @seq_parts = split /N+/ , $seq ;


		#shift(@seq_parts) ;


	#	print "Number of parts (gaps): $#seq_parts\n" ;
	
		my $start = 1 ;
		my $count = 0 ; 
		for (my $i = 0 ; $i < @seq_parts ; $i++) {
	
	#	    print "$seq_parts[$i]\n\n" ;

			my $seq_len = length($seq_parts[$i]) ;
	
			if ($seq_len < $smallest ) {
				next ;
			}
	
			#my $position_start = index($seq, $seq_parts[$i]) + 1 ;
			#my $position_end = $position_start + $seq_len - 1 ;
	
			if ( length($seq) >= $threshold ) { 
			    #print COORD "$seq_name.$i.$seq_len $position_start $position_end\n" ;
			    print OUT ">$seq_name.$count.$seq_len\n$seq_parts[$i]\n" ;
			    print READ "$seq_name.$count.$seq_len\t$seq_name\n" ;
			    $count++ ;
			}
			else {
#			    print READ2 "$seq_name.$count.$seq_len\t$seq_name\n" ;
#			    print OUT2 ">$seq_name.$count.$seq_len\n$seq_parts[$i]\n" ;
			    $count++ ;
			}


		}
	}
	else {

		
		my $seq_len = length($seq) ;

		if ( length($seq) >= $threshold) {
		    
		    #print COORD "$seq_name.noN.$seq_len 1 $seq_len\n" ;
		    print OUT ">$seq_name.noN.$seq_len\n$seq\n" ;
		    print READ "$seq_name.noN.$seq_len\t$seq_name\n" ;
		}
		elsif ( length($seq) >= $smallest ) {
		    #print OUT2 ">$seq_name.noN.$seq_len\n$seq\n" ;
 #                   print READ2 "$seq_name.noN.$seq_len\t$seq_name\n" ;

		}

	}


	
    }
	
    
    #last;


}


print "all done!\n" ; 





sub checkpath {
my $program = shift ;

	if ( `which $program` ) {
		print "found! $program path: " , `which $program` ;
	}
	else {
		print "$program path not found!\n" ;
		exit ;
	}
}

