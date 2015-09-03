#!/usr/bin/perl -w
#
#	40_walk_stage2.pl
#	Parse the ssaha output of velvet contig to reference
#
#
#


use strict;


if (@ARGV != 3) {
	print "40_walk_stage2.pl insert_size extend_bp overhang\n\n" ;
	exit ;
}


my $insert_size = shift @ARGV ;
my $extend_bp = shift @ARGV ;
my $overhang_bp = shift @ARGV ;

my $read_placed = "image.read.placed" ;
my $contig_len_file = "contigs.fa.len.txt" ;
my $insert = "image.insert_size" ;

my %contig_ins_start = ();
my %contig_ins_end = ();



# first putting the parallelised job back into the file
system("cat walks/*.walk > walk1.auto") ;
system("cat walks/*.fa > Illumina_concensus.fa") ;

# we can pipe this but let's stay this way for now.
my $stage1 = "walk1.auto" ;
my $Illumina_seq = "Illumina_concensus.fa";


# keep a track of unassembled gap
open OUT, '>', "walk2.unassembled.gaps" or die "can not create walk2.unassembled.gaps" ;
open OUT_SUM, '>', "walk2.summary" or die "can not create walk2.summary" ;

#if (-e "unassembled.gaps") {
#	open( IN, "unassembled.gaps" ) or die "Cannot open unassembled.gaps\n";
#	while (<IN>) {
#		print OUT $_ ;
#	}
#	close(IN) ;
#}







###################################################################################
#put contig length into a hash
my %contigs_len ;
my @contigs ;
my %velvet_contig ;

open( IN, "$contig_len_file" ) or die "Cannot open gzipped fastq file\n";
while (<IN>) {
	if (/(^\S+)\s+(\d+)/)  {
		$contigs_len{$1} = $2 ;
		push(@contigs, $1) ;
	}
}
close(IN);

#put illumina seq into hash
		open( IN_FASTA, "$Illumina_seq" ) or die "Cannot open velvet contigs file\n";
		my $seq_name ;
		while (<IN_FASTA>) {
			chomp ;

			if (/>(\S+)/) {
				$seq_name = $1 ;
				#print "$seq_name\n" ;
			}
			else {
				$velvet_contig{$seq_name} .= $_ ;
			}
		}
		close(IN_FASTA);

#		


###################################################################################


#use insert size
my $insert_is_file = 0 ;

if ($insert =~ /\D+/) {
	print "using the file instead!\n" ;
	open( IN, "$insert" ) or die "Cannot open gzipped fastq file\n";
	while (<IN>) {
		if (/(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/) {	
			#print "$1\t$2\t$3\t$4\n" ;
			chomp;
	
			if ( $contigs_len{$1} ) {
				if ($2 == 0 ) {
					$contig_ins_start{$1} = $insert_size ;
				}
				else {
					$contig_ins_start{$1} = $2 ;
				}

				if ($4 == 0 ) {
					$contig_ins_end{$1} = $contigs_len{$1} - $insert_size ;
				}
				else {
					$contig_ins_end{$1} = $contigs_len{$1} - $4 ;
				}
			}

			
		}
	}
	close(IN);
	$insert_is_file = 1 ;
}
if ($insert_is_file == 0 ) {
	foreach (@contigs) {
		$contig_ins_start{$_} = $insert  ;
		$contig_ins_end{$_} = $contigs_len{$_} - $insert  ;
	}
}
###################################################################################


my $possible_gaps = 0 ;
my $difficult_end = 0 ;
my $repetitive_contig = 0 ;
my $gap_closed = 0 ;

my $small_extension = 0 ;
my $big_overhang = 0 ;
my $contig_found_middle = 0 ;

my $manual_editing_needed = 0 ;


my $inversion_found = 0 ;

my %all_contigs = () ;
my %contig_extended = () ;
my %contig_closed = () ;


my $possible_gaps_check = 0 ;

my $oops = 0 ;

open( IN, "$stage1" ) or die "Cannot open walk 1 output!\n";

while (<IN>) {

	if ( /^>/ ) {
		$possible_gaps++ ;

		my @gapinfo = split /\s+/  ;
		my $contig_L = $gapinfo[1] ;
		my $contig_R = $gapinfo[3] ;

		$all_contigs{"$contig_L.$contig_R"}++ ;


		if ( $contig_L eq 'xxxxxx' ) {
			print "C:xxxxxx\tL:0\tC:$contig_R\tL:$contigs_len{$contig_R}\n" ;
		}
		elsif ( $contig_R eq 'xxxxxx' ) {
			print "C:$contig_L\tL:$contigs_len{$contig_L}\tC:xxxxxx\tL:0\n" ;
		}
		else {
			print "C:$contig_L\tL:$contigs_len{$contig_L}\tC:$contig_R\tL:$contigs_len{$contig_R}\n" ;
		}




		
		my %contigs_map = () ;
		my %contigs_solexa_coordinates = () ;
		my %contigs_reference_coordinates = () ;
		my %contigs_cigars = () ;
		my %contigs_strand = () ;


		my $small_extension_tmp = 0 ;
		my $big_overhang_tmp = 0 ;
		my $contig_middle_tmp = 0 ;

		
		while (<IN>) {
			chomp;



			
			if (/kmer/ ) {
				print "$_\n" ;
				next ;
			}
			elsif (/misc/) {
				$difficult_end++ ;

				$possible_gaps_check++ ;

				print "$_\n" ;
				print "FINISH\n\n" ;
				last; #jump out of the inside while loop
				
			}
			elsif (/cigar/) {
				#print "$_\n" ;
				
				my @cigarlines = split /\s+/ ;


	#			print "@cigarlines started!\n" ;
				
				# if the score is less than 50, then it's probably a bad match
				# next if $cigarlines[9] < 50 ;
				
				# if we have outside_range, then it's just spruious alignment
				next if /outside_range/ ;


				# exclude contigs that is matched at the contig but has extended less than the
				# specified size
				if ( $contig_L eq $cigarlines[5] ) {
					if ( $cigarlines[1] =~ /length_(\d+)/ ) {
						if ( $cigarlines[4] eq '+' && ( $1 - $cigarlines[3] ) < $extend_bp ) {
							print "misc: this match is inside contig and extends less than $extend_bp , exclude\n\n" ;
							$small_extension_tmp = 1 ;
							next ;
						}
						elsif ( $cigarlines[4] eq '-' && $cigarlines[3] < $extend_bp ) {
							print "misc: this match is inside contig and extends less than $extend_bp , exclude\n\n" ;
							#print "erm!\n";
							$small_extension_tmp = 1 ;
							next ;
						}
					}
				}
				elsif ( $contig_R eq $cigarlines[5] ) {
					if ( $cigarlines[1] =~ /length_(\d+)/ ) {
						
						if ( $cigarlines[4] eq '+' && $cigarlines[2] <  $extend_bp ) {
							print "misc: this match is inside contig and extends less than $extend_bp , exclude\n\n" ;
							$small_extension_tmp = 1 ;
							next ;
						}
						elsif ( $cigarlines[4] eq '-' && ($1 - $cigarlines[3]) < $extend_bp ) {
							print "misc: this match is inside contig and extends less than $extend_bp , exclude\n\n" ;
							$small_extension_tmp = 1 ;
							next ;
						}
					}

				}
	
		
				

		#		print "@cigarlines passed!\n";





				

				#count the number of reference contigs the velvet contig has mapped
				$contigs_map{$cigarlines[1]}{$cigarlines[5]}++ ;
				$contigs_cigars{$cigarlines[1]} .= "$_ ";
				$contigs_strand{$cigarlines[1]}{$cigarlines[4]}++ ;


			}
			else {

				#do processing here
				my $gap_close = 0 ;
				my $ambiguous = 0 ;				
				my $contig_name = '';
					

				for my $velvet_contig_tmp (keys %contigs_map) {
					#print "$velvet_contig_tmp\n" ;
					#print "number of keys:" . (scalar keys %{$contigs_map{$velvet_contig_tmp}}) . "\n" ;

					if ( (keys %{$contigs_map{$velvet_contig_tmp}}) == 2 ) {
			
						$gap_close = 1 ;
						$contig_name = $velvet_contig_tmp ;
						


						for my $contig_mapped (keys %{$contigs_map{$velvet_contig_tmp}}) {					
							if ( $contigs_map{$velvet_contig_tmp}{$contig_mapped} > 1 ) {
								print "gap closed but probably needs manual editing! \n" ;
								
								$ambiguous = 1 ;
							}
			
							

						}
						
						if ( (keys %{$contigs_strand{$velvet_contig_tmp}}) != 1 ) {
						    print "inversion found! \n" ;
						    $inversion_found++ ;
						    $ambiguous = 1 ;
						}
						

					}
					elsif ( (keys %{$contigs_map{$velvet_contig_tmp}}) > 2 ) {
						print "oops! something's wrong in...\n" ;
						exit ;
					}
					elsif ( (keys %{$contigs_map{$velvet_contig_tmp}}) == 1 ) {
						
						for ( keys %{$contigs_map{$velvet_contig_tmp}} ) {

							if ( $contigs_map{$velvet_contig_tmp}{$_} > 1 ) {
	
								print "mapped more than one place!\n" ;

								$repetitive_contig++ ;
								$ambiguous = 1 ;
							}

							

						}


					}
					
					


				}




					if ( $ambiguous == 1 ) {
						print "erm!!\n" ;
					}
					elsif ( $gap_close == 1 ) {

						$gap_closed++ ;

						print "$contig_name\t" . length($velvet_contig{$contig_name}) . "\t2\n" ;
						print "$contigs_cigars{$contig_name}\n\n" ;		

						$contig_closed{"$contig_L.$contig_R"}++ ;


					}
					else {

					    my $is_something = 0 ; 


						for my $contig_name_tmp (sort keys %contigs_cigars) { 			

							# exclude solexa contigs that are has way too long overhang...
							my @cigarlines = split /\s+/ , $contigs_cigars{$contig_name_tmp} ;
							my $contig_middle = $contigs_len{$cigarlines[5]} / 2 ;

							if ( $cigarlines[7] <= $contig_middle && $cigarlines[6] >= $overhang_bp ) {
								print "this match will replace overhang longer than $overhang_bp, exclude!\n\n" ;
								$big_overhang_tmp = 1 ;
								next ;
							}
							elsif ( $cigarlines[6] >= $contig_middle && ($contigs_len{$cigarlines[5]} - $cigarlines[7]) >= $overhang_bp ) {
								print "this match will replace overhang longer than $overhang_bp, exclude!\n\n" ;
								$big_overhang_tmp = 1 ;
								next ;
							}		

							$is_something = 1 ;
							$contig_extended{"$gapinfo[1].$gapinfo[3]"}++ ;

							print "$contig_name_tmp\t" . length($velvet_contig{$contig_name_tmp}) . "\t1\n";
							print "$contigs_cigars{$contig_name_tmp}\n\n" ;

						}




					    if ( $is_something == 0 && $big_overhang_tmp > 0 ) {
						$big_overhang++ ;
					    }
					    elsif ( $is_something == 0 && $small_extension_tmp > 0 ) {
						$small_extension++ ;
					    }
					    elsif ( $is_something == 0 && $contig_middle_tmp > 0 ) {
						$contig_found_middle++ ;
					    }
					    elsif ( $is_something == 0 ) {
						print "oops!!\n" ;
						$oops++ ;
					    }

						

					    

					    

					}
				
				

				#for (keys %contigs_solexa_coordinates) {
				#	for (keys % { $contigs_solexa_coordinates{$_} } ) {
				#		#print "$_\n" ;
				#	}
				#}


				
				$possible_gaps_check++ ;
				
				print "FINISH\n\n" ;
				last; #jump out of the inside while loop

				
				
			}
			
		}
	}








} # end of file reading


my $extend_once = 0 ;
my $extend_twice = 0 ;
my $extend_tmp = 0 ;

for (sort keys %contig_extended) {

	#print "$_\n" ;

	$extend_tmp++   if $contig_extended{$_} > 2 ;
	$extend_twice++  if $contig_extended{$_} == 2 ;
	$extend_once++  if $contig_extended{$_} == 1 ;
}


for (sort keys %all_contigs) {

	#print "$_\t" ;

	if ( $contig_extended{$_} ) { 
		#print "extended\n" ;
	}
	elsif ( $contig_closed{$_} ) { 
		#print "closed\n" ;

	}
	else {
		#print "here!\n" ;
		print OUT "$_\n" ;
	}

}
close(OUT) ;





print OUT_SUM "\n\n\n\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n" ;
print OUT_SUM "\# Stats: Total gaps: $possible_gaps\n" ;
print OUT_SUM "\# Stats: $gap_closed gap closed\n" ;
print OUT_SUM "\# Stats: $extend_once gaps extends from one side\n" ;
print OUT_SUM "\# Stats: $extend_twice gaps extends from both side\n" ;
print OUT_SUM "\# Stats: Gap produce repetitive Illumina contigs: $repetitive_contig\n" ;
print OUT_SUM "\# Stats: Gap failed to assemble at all: $difficult_end \n" ;
print OUT_SUM "\# Stats: Gap can be closed but it was inversion: $inversion_found\n" ;
print OUT_SUM "\# Stats: Gap has contigs having longer overhang than $overhang_bp : $big_overhang\n" ;
print OUT_SUM "\# Stats: Gap has contig fail within the contig ends:  $contig_found_middle\n" ;
print OUT_SUM "\# Stats: Gap has contig fail to extend more than $extend_bp bp: $small_extension\n" ;
print OUT_SUM "\# Stats: The rest is small contigs failed to extend more than $extend_bp bp at contig ends\n" ;
print OUT_SUM "\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\#\n" ;


#print "$oops\n" ;
#print "$possible_gaps_check\n" ;

