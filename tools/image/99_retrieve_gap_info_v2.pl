#!/usr/bin/perl -w
#
#	30_walk_stage1.pl
#	Local assemblies of velvet
#
#
#
#


use strict;
my $PI = `echo $$` ; chomp($PI) ;



if (@ARGV != 6) {
	print "99_retrieve_gap_info.pl old_read_placed old_seq new_read_placed new_seq out_prefix fasta_length\n" ;	

	print "example 1 (in the final directory)\n 99_retrieve_gap_info.pl ../read.placed ../contigs.fa new.read.placed new.fa final 100\n" ;

	exit ;
}

#all the parameters
my $read_placed_old  = shift @ARGV ;
my $seq_old = shift @ARGV ;
my $read_placed_new  = shift @ARGV ;
my $seq_new = shift @ARGV ;
my $out = shift @ARGV ;
my $length = shift @ARGV ;


# all the output of comparisons
open OUT, '>', "$out.stats" or die "can not create $out.stats" ;
open OUT_RAW, '>', "$out.nucmerstats" or die "can not create $out.nucmerstats" ;
open OUTFA, '>', "$out.consensus.fa" or die "can not create $out.consensus.fa" ;


my %contig_in_supercontigs_O = () ;
my %supercontigs_O = () ;
my @supercontigs_ordered_O = () ;

my %contig_in_supercontigs_N = () ;
my %supercontigs_N = () ;
my @supercontigs_ordered_N = () ;

#read in contigs in fasta format
my %contigs_fasta_O ;
my %contigs_fasta_N ;

my %contig_still_present = () ;

	# read all the files
	open( IN, "$read_placed_old" ) or die "Cannot open $read_placed_old\n";
	while (<IN>) {
	
		if (/(^\S+)\s+(\S+)/) {	
			#print "$1\t$2\n" ;

		    next if $1 eq 'xxxxxx' ;

			push @{ $contig_in_supercontigs_O{$2} } , "$1" ;
			unless ( $supercontigs_O{$2} ) {
				push(@supercontigs_ordered_O, $2) ;
			}	
			$supercontigs_O{$2}++ ;
		}
	}
	close(IN);
	
	
	open( IN, "$read_placed_new" ) or die "Cannot open $read_placed_new\n";
	while (<IN>) {
	
		if (/(^\S+)\s+(\S+)/) {	
			#print "$1\t$2\n" ;

		    next if $1 eq 'xxxxxx' ;

			push @{ $contig_in_supercontigs_N{$2} } , "$1" ;
				unless ( $supercontigs_N{$2} ) {
				push(@supercontigs_ordered_N, $2) ;
			}	
			$supercontigs_N{$2}++ ;
		    $contig_still_present{$1}++ ;
		}
	}
	close(IN);




	print "reading contigs in fasta in format....\n" ;
	open( IN, "$seq_old" ) or die "Cannot open reference fastq file: $seq_old\n";
	my $read_name = '' ;
	while (<IN>) {

	    if (/^>(\S+)/) {
		$read_name = $1 ;

		while (<IN>) {	
			if (/^>(\S+)/) {
				$read_name = $1 ;
			}
			else {
				chomp ;
				$contigs_fasta_O{$read_name} .= $_ ;
			}
		}
	    }

	}
	close(IN);



	print "reading contigs in fasta in format....\n" ;
	open( IN, "$seq_new" ) or die "Cannot open reference fastq file: $seq_new\n";
	$read_name = '' ;
	while (<IN>) {

	    if (/^>(\S+)/) {
		$read_name = $1 ;

		while (<IN>) {	
			if (/^>(\S+)/) {
				$read_name = $1 ;
			}
			else {
				chomp ;
				$contigs_fasta_N{$read_name} .= $_ ;
			}
		}
	    }

	}
	close(IN);









# stats for gaps
my $gap_closed = 0 ;
my $extended = 0 ;
my $seq_count = 0 ;
my $positive_gaps = 0 ;

my $right_ext = 0 ;
my $last_contig_name = 0 ;
my $last_supercontig_name = 0 ;

my $previous_contig_end = 0 ;
my $previous_contig_name = 0 ;






for my $supercontig ( @supercontigs_ordered_N ) {

	my @contigs_O = @{ $contig_in_supercontigs_O{$supercontig} } ;
	my @contigs_N = @{ $contig_in_supercontigs_N{$supercontig} } ;


	print "doing SC: $supercontig\n" ;
	my $count = 0 ;
	system("rm $PI.test.result") if -e "$PI.test.result" ;





	my $contig_previous = '' ;

        for ( my $i = 0 ; $i < @contigs_O ; $i++ ) {
	    print "O $contigs_O[$i]\n" ;

	    if ( $contig_still_present{$contigs_O[$i]}  ) {
		# 1 to 1 relationship
		open OUT_REF, '>', "$PI.ref.fa" or die "cannot open tmp contig file" ;
		open OUT_QUE, '>', "$PI.query.fa" or die "cannot open tmp contig file" ;
		print OUT_REF ">$contigs_O[$i]\n$contigs_fasta_N{$contigs_O[$i]}\n" ;
		print OUT_QUE ">$contigs_O[$i]\n$contigs_fasta_O{$contigs_O[$i]}\n" ;
		close(OUT_REF) ;
		close(OUT_QUE) ;

		$contig_previous = "$contigs_O[$i]" ;

		system("nucmer --prefix=$PI.ref_qry $PI.ref.fa $PI.query.fa") ;
		system("show-coords -rclH $PI.ref_qry.delta | awk '\$16 > 70' >> $PI.test.result") ; 

	    }
	    else {

		open OUT_REF, '>', "$PI.ref.fa" or die "cannot open tmp contig file" ;
                open OUT_QUE, '>', "$PI.query.fa" or die "cannot open tmp contig file" ;
		print OUT_REF ">$contig_previous\n$contigs_fasta_N{$contig_previous}\n" ;
                print OUT_QUE ">$contigs_O[$i]\n$contigs_fasta_O{$contigs_O[$i]}\n" ;
		close(OUT_REF) ;
                close(OUT_QUE) ;

		print "$contigs_O[$i] link to $contig_previous\n" ;

                system("nucmer --prefix=$PI.ref_qry $PI.ref.fa $PI.query.fa") ;
		system("show-coords -rclH $PI.ref_qry.delta | awk '\$16 > 70' >> $PI.test.result") ;

	    }


        }







	open (RESULT_RAW, "$PI.test.result") or die "can not pipe result of nucmer..\n" ;
	

	my %old_contig_present = () ;
	my $contig_order = 0 ;



	while (<RESULT_RAW>) {
		chomp ;
		print OUT_RAW "$_\n" ;
		print "$_\n" ;
		

		my @line = split /\s+/ , $_ ;
		shift(@line) ;

		$contig_order = 0 if $line[17] eq $line[18] ;

		#print OUT "0: $line[0], 9: $line[9]\n" ;

    		#[S1] [E1] | [S2] [E2] |  LEN 1] [LEN 2] | [% IDY]  | [LEN R] [LEN Q]  | [COV R] [COV Q] | [TAGS]
		# 0    1   2  3    4   5   6       7     8    9     10  11      12     13  14      15    16  17  18  


		
		if ( $contig_order == 0 ) {




			if ( $right_ext != 0 ) {

                        # print out 3ext fasta
			    my $left_pos = $previous_contig_end  -1 - $length ;
			    my $fasta_seq = substr ( $contigs_fasta_N{"$previous_contig_name"} , $left_pos  ) ;


#			      print "$previous_contig_name\n" ;

			    print OUTFA ">$previous_contig_name.3ext.from" . ($left_pos + 1 ) . "\n" ;
			    print OUTFA "$fasta_seq\n" ;

			    print OUT "$supercontig\textend\t3end\t$right_ext $previous_contig_name.3ext.from" . ($left_pos + 1 ) . "\n" ;

			    $seq_count++ ;
			    $extended++ ;

			}



			if ( ( $line[0] - 1 ) != 0 ) {

				my $left_pos = 0 ;
                                my $parsed_seq_len = $line[0] - 1 + $length ;
				my $fasta_seq = substr ( $contigs_fasta_N{"$line[17]"} , $left_pos , $parsed_seq_len ) ;

				print OUTFA ">$line[17].5ext.$parsed_seq_len\n" ;
				print OUTFA "$fasta_seq\n" ;

				print OUT "$supercontig\textend\t5end\t". ( $line[0] - 1 )   .  " $line[17].5ext.$parsed_seq_len\n" ;

				$seq_count++ ;
				$extended++ ;

			}


		}


		my $l_overhang = $line[3] - 1 ;
		my $r_overhang = $line[12] - $line[4] ;

		if ($contig_order != 0 ) {

			my $gap_len = $line[0] - $previous_contig_end -1 ;

			if ( $gap_len > 0 ) {


				my $left_pos = $previous_contig_end -1 - $length ;
				my $parsed_seq_len = $gap_len  +  $length + $length ;

                                print OUT "$supercontig\tMerge\t+ve\t$gap_len $previous_contig_name.$line[18].$parsed_seq_len\n" ;

				my $fasta_seq = substr ( $contigs_fasta_N{"$line[17]"} , $left_pos , $parsed_seq_len ) ;

#				print "$left_pos, $parsed_seq_len\n" ;


				print OUTFA ">$previous_contig_name.$line[18].$parsed_seq_len\n" ;
				print OUTFA "$fasta_seq\n" ;

				$positive_gaps++ ;
				$gap_closed++ ;
				$seq_count++ ;
			}
			else {


				print OUT "$supercontig\tMerge\t-ve\t$gap_len\n" ;

				$gap_closed++ ;

			}

			

		}

		print OUT "$supercontig\tcontig\t$line[17]\t$line[11]\t$line[18]\t$line[12]\t$line[14]\t$line[0]\t$line[1]\t$l_overhang\t$r_overhang\n" ;

		
		$right_ext = $line[11] - $line[1] ;
		$contig_order++ ;
		$previous_contig_name = $line[17] ;
		$previous_contig_end = $line[1] ;
		$last_contig_name = $line[18] ; 
		$last_supercontig_name = $supercontig ;

	}


	#last;






}


if ( $right_ext != 0 ) {

                        # print out 3ext fasta
    my $left_pos = $previous_contig_end  -1 - $length ;
    my $fasta_seq = substr ( $contigs_fasta_N{"$previous_contig_name"} , $left_pos  ) ;

    print OUTFA ">$last_contig_name.3ext.from" . ($left_pos + 1 ) . "\n" ;
    print OUTFA "$fasta_seq\n" ;

    print OUT "$last_supercontig_name\textend\t3end\t$right_ext $last_contig_name.3ext.from" . ($left_pos + 1 ) . "\n" ;

    $seq_count++ ;
    $extended++ ;

}






print "\n\n\n" ;
print "Summary of comparisons between final consensus and original contigs.fa\n" ;
print "Number of gaps closed: $gap_closed\n" ;
print "Number of +ve gaps:    $positive_gaps\n" ;
print "Number of extensions:  $extended\n" ;
print "A total of $seq_count sequences presented in $out.consensus.fasta\n" ;


print OUT "\n\n\nSummary of comparisons between final consensus and original contigs.fa\n" ;
print OUT "Number of gaps closed: $gap_closed\n" ;
print OUT "Number of +ve gaps:    $positive_gaps\n" ;
print OUT "Number of extensions:  $extended\n" ;
print OUT "A total of $seq_count sequences presented in $out.consensus.fasta\n" ;


system("rm $PI.ref.fa $PI.query.fa $PI.ref_qry.*") ;

close(OUT) ;
close(OUTFA) ;






 
