#!/usr/bin/perl
#
#	generating a concensus sequence after deciding which solexa contigs to merge
#	walk_stage3 stand alone version
#
#	author: jit
#	last updated 30th October
#
#	function added: some additional codes of taking care of smaller contigs.
#
#
use strict;



#if (@ARGV != 1) {
#	print "$0 scaff.num \n\n" ;
#	exit ;
#}




#my $scaffnum = shift @ARGV ;


my $stage2 = 'walk2.auto';
my $output = 'new';






my $contig_fa = "image.contigs.fa";
my $read_placed = "image.read.placed";


my $illumina_concensus = 'Illumina_concensus.fa' ;

#opening output of all sorts
my $out_fa ; my $out_insert ; my $out_placed ; my $out_report ;
my $out_fa_small ; 

open $out_fa, ">" , "$output.fa" or die "can not open $output.fa\n" ;
open $out_fa_small , ">", "$output.small.fa" or die "cannot open $output.small.fa\n" ; 


open $out_insert, ">" , "$output.insert_size" or die "can not open $output.stats\n" ;
open $out_placed, ">" , "$output.read.placed" or die "can not open $output.read.placed\n" ;
open $out_report, ">" , "$output.report" or die "can not open $output.read.placed\n" ;


#parameters
my %contig_seq = () ;                                                                                                                                                                 
# flags for scaffolding...
my %contig_start_pos ;
my %contig_start_pos_solcontig ;
my %contig_end_pos ;
my %contig_pointers ;
my @contig_ordered ;
my %contigs_merged_into ;

my %contigs_3end_extended = () ;
my %contigs_5end_extended = () ;


print "1. reading original contig.fa file!..." ;
open (IN, "$contig_fa") or die "oops! can not open $contig_fa\n" ;
my $seq_name = $1 ;
my $seq = '' ;


while (<IN>) {
    if (/^>(\S+)/) {
	$seq_name = $1 ;

	while (<IN>) {
	    chomp ;

	    if (/^>(\S+)/) {
		#print "$seq_name\n" ;

		$contig_seq{$seq_name} = $seq ;
		$seq_name = $1 ;
		$seq = '' ;
	    }
	    else {
		$seq .= $_ ;
	    }

	}

    }
}

$contig_seq{$seq_name} = $seq ;
close(IN) ;




#used to process the final output
for ( keys %contig_seq ) {
	$contig_start_pos{$_} = 0 ;
	$contig_start_pos_solcontig{$_} = 0 ;
	$contig_end_pos{$_} = 0 ;
	$contig_pointers{$_} = $_ ;
	push @{ $contigs_merged_into{ $_ } } , $_ ;
}


# parse the velvet contig
  my %velcontig_seq = () ;
  open my $VEL, "< $illumina_concensus" or die "can not open $illumina_concensus!\n" ;
  while (<$VEL>) {
	if (/^>(\S+)/) {
		    my $seq_name = $1 ;
		    my $seq = <$VEL> ;
		    chomp($seq) ;
		    $velcontig_seq{$seq_name} = $seq ;
	}
  }
  close($VEL) ;


#start parsing the output
open( IN, "$stage2" ) or die "Cannot open$stage2\n";

#debug purposes
my $tmp_count = 0 ;








while (<IN>) {

	# if a gap is found
	if (/^C:(\S+)\s+L:(\d+)\s+C:(\S+)\s+L:(\d+)/) {

		print "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n" ;
		print "\nSTART:reading $1 and $3..\n" ;
		
		my $contigs = "$1.$3" ;
		my $contig_L = $1; 
		my $contig_R = $3; 

		#debug purpose
		$tmp_count++ ;
	
		#to jump out of particular contigs
		while (<IN>) {
		       
			#start of the contig
			if (/($contigs.NODE_\d+_length_)(\S+)\s+(\d+)\s+(\d+)/) {

			  


	
				# some stats
				my $vel_contig = "$1$2" ;
				my $vel_contig_length = length($velcontig_seq{$vel_contig}) ;
				my $ismerged = $4 ;
								


				# if the contigs are joineds
				if ( $ismerged == 2) {
						
					
					print "MERGING 2 contigs: $contig_L, $contig_R \n" ;

					my @cigarline = split /cigar/, <IN> ;
					shift(@cigarline) ; #remove the empty line
				
					my $l_pos_start ;
					my $r_pos_start ;
					my $l_pos_end ;
					my $r_pos_end ;
					my $l_strand ;
					my $r_strand ;

					my $l_ref_start ;
					my $l_ref_end ;
					my $r_ref_start ;
					my $r_ref_end ;

					my @l_matched ;
					my @r_matched ;

					my $l_found = 0 ;
					my $r_found = 0 ;

					print "\t2 cigar files:\n" ;
					foreach(@cigarline) {
						my @cigar_result = split /\s+/ , $_ ;
						print "\t@cigar_result\t" ;
						print @cigar_result . "\n" ;
	
						if ($cigar_result[5] eq "$contig_L") {

							$l_strand = $cigar_result[4] ;
							$l_pos_start = $cigar_result[2] ;
							$l_pos_end = $cigar_result[3] ;
							$l_ref_start = $cigar_result[6] ;
							$l_ref_end = $cigar_result[7] ;
							@l_matched = @cigar_result[10..$#cigar_result] ;
							$l_found = 1 ;
						}
						elsif ( $cigar_result[5] eq "$contig_R") {
							
							$r_strand = $cigar_result[4] ;
							$r_pos_start = $cigar_result[2] ;
							$r_pos_end = $cigar_result[3] ;
							$r_ref_start = $cigar_result[6] ;
							$r_ref_end = $cigar_result[7] ;
							@r_matched = @cigar_result[10..$#cigar_result] ;
							$r_found = 1 ;
						}
						else {
							print "something wrong here!!! hehe?\n" ; exit ;
						}

					}
					print "\n" ;

			

					if ($l_strand eq '-') {

					    $velcontig_seq{$vel_contig} = reverse($velcontig_seq{$vel_contig}) ;

						my $reversecomplemented = $velcontig_seq{$vel_contig} ;
						$reversecomplemented =~ tr/ATCGatcg/TAGCTAGC/ ; 				
						$velcontig_seq{$vel_contig} = $reversecomplemented ;
							 				
						$l_pos_start= $vel_contig_length - $l_pos_start +1;
						$r_pos_start= $vel_contig_length - $r_pos_start +1;
						$l_pos_end= $vel_contig_length - $l_pos_end +1 ;
						$r_pos_end= $vel_contig_length - $r_pos_end +1 ;		
					}


					

					
					

					# generate consensus ; insert indels here (i.e., add gaps)
					print "generate consensus!....\n" ;

					print "l_pos_start = $l_pos_start\n" ;
					print "l_pos_end = $l_pos_end\n" ;
					print "r_pos_start = $r_pos_start\n" ;
					print "r_pos_end = $r_pos_end\n" ;

					print "l_ref_start = $l_ref_start\n" ;
					print "l_ref_end = $l_ref_end\n" ;
					print "r_ref_start = $r_ref_start\n" ;
					print "r_ref_end = $r_ref_end\n" ;


					#where the velvet contig should end
					my $contig_to_insert = substr($velcontig_seq{$vel_contig}, ($l_pos_start-1),  ( $r_pos_end - $l_pos_start +1) ) ;

					my $left_contig_length = length( substr($contig_seq{$contig_pointers{$contig_L}}, 0,($l_ref_start-1-1 + $contig_start_pos{  $contig_L  }      )) ) ;


						$contig_seq{ $contig_pointers{$contig_L} } = substr($contig_seq{$contig_pointers{$contig_L}}, 0,($l_ref_start-1 + $contig_start_pos{  $contig_L  }      )) 
							. $contig_to_insert 
							. substr($contig_seq{$contig_R}, $r_ref_end) ;


					if ( $l_pos_end > $r_pos_start ) {
						print "\n\noverlap found!\n" ;




						print $out_report "MERGED: $contig_pointers{$contig_L} <- $contig_R\t" ;
						print $out_report "Solexa contig position: overlap was found from $l_ref_start of $contig_L to $r_ref_end of $contig_R \n" ;

						$contig_start_pos{$contig_R} =  $left_contig_length + $r_pos_start - ($r_ref_start - 1);
						print "new contig_R start pos: " . $contig_start_pos{$contig_R} . "\n" ;

					}
					else {


						print $out_report "MERGED: $contig_pointers{$contig_L} <- $contig_R\n" ;		

						#debugging		
						my $sol_contig_start = $l_ref_start-1 + $contig_start_pos{  $contig_L  } ;
						$contig_start_pos{$contig_R} = length($contig_seq{ $contig_pointers{$contig_L} }) - length(substr($contig_seq{$contig_R}, $r_ref_end)) - $r_ref_end ;

						print "new contig_R start pos: " . $contig_start_pos{$contig_R} . "\n" ;
						
		

					}

					

						$contig_pointers{$contig_R} =  $contig_pointers{$contig_L} ;
						push @{ $contigs_merged_into{ $contig_pointers{$contig_L} } } , $contig_R ;
	


					
						print "MERGING DONE!\n" ; last;
					
				}
				# if there's only one match
				else {

#my %contig_start_pos ;
#my %contig_pointers ;

					print "EXTENSION found!...\n" ;

					my $ismerged = 0 ;

					my $pos_start ;
					my $pos_end ;
					my $strand ;
					my $ref_start ;
					my $ref_end ;
					my $l_found = 0 ;
					my $r_found = 0 ;

					my $cigar_result_tmp = <IN> ;
					chomp($cigar_result_tmp) ;
					my @cigar_result = split /\s+/ , $cigar_result_tmp ;

					print "cigar line: @cigar_result\t" ;
					print @cigar_result . "\n" ;
	
					print "contigL: $contig_L contigR: $contig_R\n" ; 

					if ($cigar_result[5] eq "$contig_L") {
						print "EXTENSION found in $contig_L\n" ;
						$l_found = 1 ;
						if ( $contig_pointers{$contig_L} ne $contig_L  ) {
							print "it has been merged previously!\n" ;
							$ismerged = 1 ;
						}

					}
					else {
						print "EXTENSION found in contig_R: $contig_R\n\n" ;
						$r_found = 1 ;
					}

	
					$strand = $cigar_result[4] ;
					$pos_start = $cigar_result[2] ;
					$pos_end = $cigar_result[3] ;
					$ref_start = $cigar_result[6] ;
					$ref_end = $cigar_result[7] ;
							

					#reverse the sequence
					if ($strand eq '-') {
						$velcontig_seq{$vel_contig} = reverse( $velcontig_seq{$vel_contig}) ;
						my $reversecomplemented = $velcontig_seq{$vel_contig};
						print "before: $reversecomplemented\n" ; 

						$reversecomplemented =~ tr/ATCGatcg/TAGCTAGC/ ;


						$velcontig_seq{$vel_contig} = $reversecomplemented ;
							 
						print "old position start and end: $pos_start $pos_end\n" ;					
						$pos_start= $vel_contig_length - $pos_start +1;
						$pos_end= $vel_contig_length - $pos_end +1 ;
						print "new position start and end: $pos_start $pos_end\n" ;					


					}


					# generate consensus ; insert indels here (i.e., add gaps)
					print "generate consensus!...\n" ;

					if ( $l_found == 1 ) {
						if ( $ismerged == 0 ) {

							if ( $ref_end == length($contig_seq{$contig_L}) && $ref_start == 1 ) {
								print "whole contig has been covered and extended!!!\n" ;
								$contig_seq{$contig_L} = $velcontig_seq{$vel_contig}  ;								
								$contig_end_pos{$contig_L} =  length($velcontig_seq{$vel_contig});

							}
							elsif ( (length($contig_seq{$contig_L}) - $contig_start_pos{$contig_L} - $ref_end ) < ($ref_start - 1) ) {
								print "(extend on the 3' side...) ref_start: $ref_start\t contig_start_pos: $contig_start_pos{$contig_L}\n" ;
								print $out_report "EXTEND: $contig_L \t 3' end from " . (  $ref_start + $contig_start_pos{$contig_L} - $pos_start ) . "\n";
								$contig_seq{$contig_L} = substr($contig_seq{$contig_L}, 0,(  $ref_start + $contig_start_pos{$contig_L} - $pos_start )) . substr($velcontig_seq{$vel_contig},0)  ;								
								$contig_end_pos{$contig_L} =  length($velcontig_seq{$vel_contig});


								if ( $contigs_3end_extended{$contig_L} ) {
									print "already extended!! next\n" ;
									next ;
								}


							}
							else {
								print "merging on 5 side!??!, do something here.....\n" ;

								if ( $contigs_5end_extended{$contig_L} ) {
									print "already extended!! next\n" ;
									next ;
								}

								print "(extend on the 5' side...) ref_start: $ref_start\t contig_start_pos: $contig_start_pos{$contig_L}\n" ;
								$contig_seq{$contig_L} = substr($velcontig_seq{$vel_contig},0,$pos_end) . substr($contig_seq{$contig_L}, $ref_end) ;				
								$contig_start_pos{$contig_L} =  length(substr($velcontig_seq{$vel_contig},0,$pos_end)) - $ref_end  ;
								$contig_start_pos_solcontig{$contig_L} = length($velcontig_seq{$vel_contig}) ;					
			


							}

						}
						else {											

								print "(extend on the 3' side...) ref_start: $ref_start\t contig_start_pos: $contig_start_pos{$contig_L}\n" ;
								print $out_report "EXTEND: $contig_pointers{$contig_L} \t 3' end from " . (  $ref_start + $contig_start_pos{ $contig_L } - $pos_start ) . "\n";
								$contig_seq{ $contig_pointers{$contig_L} } = substr($contig_seq{ $contig_pointers{$contig_L} }, 0,(  $ref_start + $contig_start_pos{ $contig_L } - $pos_start )) . substr($velcontig_seq{$vel_contig},0)  ;								
								$contig_end_pos{$contig_L} =  length($velcontig_seq{$vel_contig});							
							

						}

					}
					elsif ( $r_found == 1 ) {
					    print $out_report "EXTEND: $contig_R \t 5' end to $ref_end\n";						


	
						if ( $ref_end == length($contig_seq{$contig_R}) && $ref_start == 1 ) {
								print "whole contig has been covered and extended!!!\n" ;
								$contig_seq{$contig_R} = $velcontig_seq{$vel_contig}  ;								
								$contig_end_pos{$contig_R} =  length($velcontig_seq{$vel_contig});

						}
						# if it's a small contig and actually extends on the 3 end...
						elsif ( (length($contig_seq{$contig_R}) - $ref_end ) < ($ref_start - 1) ) {
							print "(extend on the 3' side!!! at 3' contig, prob. small contigs ) ref_start: $ref_start\t contig_start_pos: $contig_start_pos{$contig_R}\n" ;
							$contig_seq{$contig_R} = substr($contig_seq{$contig_R}, 0, $ref_start) . substr($velcontig_seq{$vel_contig},0)   ;				

							$contigs_3end_extended{$contig_R}++ ;

	
						}
						else {
							print "(extend on the 5' side...) ref_start: $ref_start\t contig_start_pos: $contig_start_pos{$contig_R}\n" ;
							$contig_seq{$contig_R} = substr($velcontig_seq{$vel_contig},0,$pos_end) . substr($contig_seq{$contig_R}, $ref_end) ;				
							$contig_start_pos{$contig_R} =  length(substr($velcontig_seq{$vel_contig},0,$pos_end)) - $ref_end  ;
							$contig_start_pos_solcontig{$contig_R} = length($velcontig_seq{$vel_contig}) ;	

							$contigs_5end_extended{$contig_R}++ ;
						}




					}								
					print "generate consensus!...done!\n\n" ;




					print "EXTENSION done!\n" ;


					#print "" . $contig_seq{'00015.scaffold23.1.size518936.14.48778'} . "\n" ; 


				}	# end of either join or just extend loop


				
			}
			elsif (/FINISH/) {
				 print "\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n\n\n" ;
				 last ;
			}


		}
	    }
	
	#last if $tmp_count == 100;

}
close(IN) ;



##################################################################################################

print "2: NOW PRODUCE THE FINAL AFG FILE WITH JOINED CONTIGS!\n" ;
print "output file: $output.fa , $output.insert_size $output.read.placed \n\n" ;

# output all merged contigs.
my $test_contig = 'contig_855' ;
my $test_mode = 0 ;

print "\n\n" ;






my %contig_too_small_so_excluded = () ; 


# print out all the contig sequences
for my $merged_contigs (keys %contig_pointers) {

	next if $merged_contigs ne $contig_pointers{$merged_contigs} ;

	$merged_contigs = $test_contig if $test_mode == 1;
	
	my $merged_contigs_tmp = $merged_contigs ;

	#stats and new concensus

	if ( length($contig_seq{$merged_contigs}) >= 300 ) {

		print $out_fa ">$merged_contigs\n" or die "truncated printing!\n";
		print $out_fa $contig_seq{$merged_contigs} . "\n" or die "truncated printing!\n" ;

		print ">$merged_contigs\n\n" ; 
		print $contig_seq{$merged_contigs} . "\n\n" ; 

		#print $out_report "$merged_contigs contains: @{ $contigs_merged_into{$merged_contigs} } \n" ;
		my $end_contig = @{ $contigs_merged_into{$merged_contigs} }[$#{ $contigs_merged_into{$merged_contigs} }] ;
		print $out_insert "$merged_contigs\t$contig_start_pos_solcontig{$merged_contigs}\t$end_contig\t$contig_end_pos{$end_contig}\n" ;


#		last ; 
	}
	else {
	    $contig_too_small_so_excluded{$merged_contigs}++ ; 

	    print $out_fa_small ">$merged_contigs\n" ; 
	    print $out_fa_small $contig_seq{$merged_contigs} . "\n" or die "truncated printing!\n" ;
	}

	#$iid_for_contig++ ;
	last if $test_mode == 1;
	

}




################################################################################################################
#final output 

my %contigs_in_supercontigs ;
my %contigs_found ;
my %supercontigs_found ;
my @supercontigs_ordered ;

print "now, creating new read.placed!...\n" ;
open my $old_read_placed, "< $read_placed" or die "can not open $read_placed!\n" ;
print "readplaced file is : $read_placed\n" ;

while(<$old_read_placed>) {
    my @read_placed_split = split /\s+/ , $_ ;

    if ( $contigs_found{$read_placed_split[0]}  ) {

    }
    else {
	$contigs_in_supercontigs{$read_placed_split[1]} .= "$read_placed_split[0] " ;
	$contigs_found{$read_placed_split[0]}++ ;
	if ( $supercontigs_found{$read_placed_split[1]} ) {

	}
        else {
	   push(@supercontigs_ordered, $read_placed_split[1]) ;
	   $supercontigs_found{$read_placed_split[1]}++ ;
	}
	 
    }
}
close($old_read_placed) ;


my $scaffold_count = 0 ;


	foreach my $supercontig (@supercontigs_ordered) {

	  

	   my @contigs = split /\s+/ ,  $contigs_in_supercontigs{$supercontig};

	    #print "doing SC: $supercontig\n" ;


		for ( my $i = 0 ; $i < ( $#contigs + 1) ; $i++ ) {

			if ($i == 0 ) {
				print $out_placed "xxxxxx\t $supercontig\n" ;
			       
			}


			if ($i == $#contigs ) {	    
				if ( $contigs[$i] eq $contig_pointers{$contigs[$i]} ) {

				    if ( $contig_too_small_so_excluded{$contigs[$i]} ) {
                        print "$contigs[$i] too small!\n" ;
				    }
				    else {
                        print $out_placed "$contigs[$i]\t $supercontig\n" ;
					    print $out_placed "xxxxxx\t $supercontig\n" ;
					
				    }
				}
				else {
					print $out_placed "xxxxxx\t $supercontig\n" ;
				}
			
				next ;
			}

			next if $contigs[$i] ne $contig_pointers{$contigs[$i]} ;

			if ( $contig_too_small_so_excluded{$contigs[$i]} ) {
			    
			    print "$contigs[$i] too small!\n" ;
			
			}
		        else {
			    print $out_placed "$contigs[$i]\t$supercontig\n" ;
			}
			
		}


	   $scaffold_count++ ;


	}


close($out_placed) ;

#system("cp $output.read.placed $output.read.placed.backup") ;
#system("01_prepare_new_read_placed.pl $output.read.placed") ;


print "all done!\n\n" ;


















