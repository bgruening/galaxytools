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
my $path ='' ; 

if ($0 =~ /(\S+)\/\w+\.pl/ ) {
    $path = $1 ;
}

my $velvetg = "$path/velvetg ";
my $velveth = "$path/velveth ";
my $ssaha = "$path/ssaha2 " ;
my $checkvelvet = "$path/check_velvet_contig.pl " ;

my $random_no = int(rand(1000000)) ;



if (@ARGV != 5) {
	print "$0 jobarray_id insert_size toexclude overhang vel_ins_len \n" ;	
	exit ;
}

#all the parameters
my $job_id = shift @ARGV ;

# sorting out different versions
my $read_placed = "read_placeds/$job_id" ;
my @kmers = qw / 27 31 41 51 / ;

if ($job_id =~ /kmer(\d+)/) { # version 1
	$read_placed = 'image.read.placed' ;
	$job_id = 1 ;
	my $kmer = $1 ;
	@kmers = ("$kmer") ;
}
my $insert_size = shift @ARGV ;

# make a few more directories
system("mkdir tmp/$PI") ;

#velvet options

my $toexclude = shift @ARGV ;
my $overhang = shift @ARGV ;

my $velvet_insert_size = shift @ARGV ;
my $file_prefix = sprintf("%06s", $job_id) ;


# walk1 output and concensus sequences
open OUT, '>', "walks/$file_prefix.walk" or die "can not create walk1 output" ;
open OUTFA, '>', "walks/$file_prefix.fa" or die "can not create Illumina_concensus.fa " ;





my $insert_len = "image.insert_size" ;
my $contigs_seq = "image.contigs.fa" ;
my $contig_len_file = "contigs.fa.len.txt" ;
my $file_location = "bridge" ;
#my $velvet_location = "velvet_$kmer_size.auto" ;




###################################################################################
#put contig length into a hash
my $insert_is_file = 0 ;
my %contig_ins_start = ();
my %contig_ins_end = ();
my @contigs_ordered ;
my %contigs_len ;

# for the singlet contig

my %Illumina_contigs = () ;

#read in the contig file
open( IN, "$contig_len_file" ) or die "Cannot open $contig_len_file\n";
while (<IN>) {
	if (/(^\S+)\s+(\d+)/) {
		push (@contigs_ordered, $1) ;
		$contigs_len{$1} = $2 ;
	}
}
close(IN);

$contigs_len{'xxxxxx'} = 'xxxxxx' ;


###################################################################################



#read in contigs in fasta format
my %contigs_fasta ;

	print "reading contigs in fasta in format....\n" ;
	open( IN, "$contigs_seq" ) or die "Cannot open reference fastq file: $contigs_seq\n";
	my $read_name = '' ;
	my $read_seq = '' ;

	while (<IN>) {
	    if (/^>(\S+)/) {
		$read_name = $1 ;
		$read_seq = ">$1\n" ;
		
		while (<IN>) {

			if (/^>(\S+)/) {
				$contigs_fasta{$read_name} = $read_seq ;     

				$read_name = $1 ;
				$read_seq = ">$1\n" ;

			}
			else {
				$read_seq .= $_ ;
			}


		}

	    }
	}
	$contigs_fasta{$read_name} = $read_seq ;    

	close(IN);

###################################################################################

# record the contig start and end position
if ($insert_len =~ /\D+/) {
	print "reading file: $insert_len\n" ;
	open( IN, "$insert_len" ) or die "Cannot open $insert_len\n";
	while (<IN>) {
		if (/(\S+)\s+(\d+)\s+(\S+)\s+(\d+)/) {	
			#print "$1\t$2\t$3\t$4\n" ;
	
			chomp;
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
	close(IN);
	$insert_is_file = 1 ;
}

###################################################################################


# a hash of unassembled contigs
my %unclosable_contigs = () ;
print "read unassembled gap in previous iteration...\n" ;
if (-e "unassembled.gaps") {
	open( IN, "unassembled.gaps" ) or die "Cannot open unassembled.gaps\n";
	while (<IN>) {
		chomp ;
		$unclosable_contigs{$_}++ ;
	}
	close(IN) ;
}



####################################################################################


# sort the contig and supercontig
open( IN, "$read_placed" ) or die "Cannot open $read_placed\n";

my %contig_in_supercontigs = () ;
my %supercontigs = () ;
my @supercontigs_ordered = () ;

my %contig_present = () ;

while (<IN>) {

    print "hmm $_" ;

	if (/(^\S+)\s+(\S+)/) {	
		#print "$1\t$2\n" ;
		push @{ $contig_in_supercontigs{$2} } , "$1" ;
		unless ( $supercontigs{$2} ) {
			push(@supercontigs_ordered, $2) ;
		}	
		$supercontigs{$2}++ ;
		$contig_present{$1}++ ;
	}
}



####################################################################################


# read the fasta of contig ends into hash
print "read illumina fasta at ends into hash..\n\n" ;
open( IN, "final.list" ) or die "Cannot open the list of fastas at contig ends\n";

my %reads_in_ends = () ;

while (<IN>) {



    	my @read_in_end = split /\s+/ ;

	unless ( $contig_present{$read_in_end[0]} ) {
		next ;
	}


    #output is like this
    # contig_0.1.30000        Lend    Read.Salmonella_paratyphi_A_chromosome.276.47   ggtaacggtgcgggctgacgcgtacaggaaacacagaaaaaagcccgcacctgaacagtgcgggcttttttttcga    aattaccaaccatctggtggcgatgattgaaaaaactatcggcggtcaggatgctttgccgaatatcagcgatgcc


    $reads_in_ends{$read_in_end[0]}{$read_in_end[1]}{$read_in_end[2]}{"1"} = $read_in_end[3] ;
    $reads_in_ends{$read_in_end[0]}{$read_in_end[1]}{$read_in_end[2]}{"2"} = $read_in_end[4] ;


}
close(IN) ;



####################################################################################









print "starting to walk...\n" ;



#my $tmp_SC = 0 ;


for my $supercontig ( @supercontigs_ordered ) {

my @contigs = @{ $contig_in_supercontigs{$supercontig} } ;

print "doing SC: $supercontig\n" ;

my $count = 0 ;

# for every contig pairs in supercontig
for ( my $i = 0 ; $i < $#contigs ; $i++ ) {

	my $contig_L = $contigs[$i] ;
	my $contig_R = $contigs[$i+1] ;
	my $contig_L_len = $contigs_len{$contig_L} ;
	my $contig_R_len = $contigs_len{$contig_R} ;

	print "doing $contig_L\t$contigs_len{$contig_L}\t$contig_R\t$contigs_len{$contig_R}...\n" ;
	print OUT ">\t$contig_L\t$contigs_len{$contig_L}\t$contig_R\t$contigs_len{$contig_R}\n" ;


	#exclude contigs if it's unable to assemble in previous iteration
	if ( $unclosable_contigs{"$contig_L.$contig_R"} ) {
		print OUT "misc: unclosable in previous iteration\n\n" ;
		next ;
	}


	#velvet stage

	# check 1: too few or too many reads
	my $read_mapped_at_contig_end =  ( keys (%{$reads_in_ends{"$contig_L"}{"Rend"}}) + keys (%{$reads_in_ends{"$contig_R"}{"Lend"}}) ) * 2;
	if ( $read_mapped_at_contig_end < 10 || $read_mapped_at_contig_end > 200000 ) {
			print "Too few or too many reads: $read_mapped_at_contig_end, skip..\n" ;
			print OUT "misc: Zero contigs were produced! (too few or too many reads: $read_mapped_at_contig_end)\n\n" ;
			next ;
	}

	
	#velveth stage
	foreach my $kmer (@kmers) {
	
		# pipe the fasta into velveth	
		if ( open(VELH, "| $velveth velvet_$kmer.auto/$contig_L.$contig_R $kmer -fasta -shortPaired -") ) {

			for my $read ( keys %{$reads_in_ends{"$contig_L"}{"Rend"}}  ) {
					print VELH ">$read/1\n" ;
					print VELH $reads_in_ends{"$contig_L"}{"Rend"}{$read}{"1"} . "\n";
					print VELH ">$read/2\n" ;
					print VELH $reads_in_ends{"$contig_L"}{"Rend"}{$read}{"2"} . "\n";
		
			}
		
			for my $read ( keys %{$reads_in_ends{"$contig_R"}{"Lend"}}  ) {
					print VELH ">$read/1\n" ;
					print VELH $reads_in_ends{"$contig_R"}{"Lend"}{$read}{"1"} . "\n";
					print VELH ">$read/2\n" ;
					print VELH $reads_in_ends{"$contig_R"}{"Lend"}{$read}{"2"} . "\n";
		
			}
		
			close(VELH) || die "can't close velveth: $!";

		}
		else {

			print "velveth did not open for some reason, skip this\n" ;
			print OUT "kmer$kmer: velveth failed to run\n\n" ;
		}

	}


	my %velvet_contig_result = () ;
	my $num_reads_used = 0 ;
	my $maxcontig_len = 0 ;
	my $kmerchoice = 31 ;


	# velvetg stage
	foreach my $kmer (@kmers) {

		# velvetg stage
		# system("rm tmp/$PI/$PI.$random_no.velg.out") if -e "tmp/$PI/$PI.$random_no.velg.out" ;

		# pipe it
		if ( -d "velvet_$kmer.auto/$contig_L.$contig_R" ) {
			system("$velvetg velvet_$kmer.auto/$contig_L.$contig_R -scaffolding no -very_clean yes -exp_cov auto -min_contig_lgth 250 -ins_length $velvet_insert_size ") ;
			
		}
		else {
			print OUT "kmer$kmer: velveth failed to run, so skip velvetg\n\n" ;
			next ;
		}	
	
		open (VELG, "velvet_$kmer.auto/$contig_L.$contig_R/Log") or die "can not open file: $!\n" ;
		# parse out result of velvet
		# pipe it straight into memory
		while (<VELG>) {
			chomp ;

			if (/Final graph/) {
				print OUT "kmer$kmer: $_\n" ;

				if (/max (\d+).+using (\d+)\// ) {
					my $tmp_max_contig_len = $1 + $kmer - 1;
					my $tmp_num_reads_mapped = $2 ;
	
					if ( $tmp_max_contig_len > $maxcontig_len ) {
							$num_reads_used = $tmp_num_reads_mapped ;
							$kmerchoice = $kmer ;
							$maxcontig_len = $tmp_max_contig_len ;
					}
				}
			}
			#elsif ( /Estimated Coverage/ ) {
			#	print OUT "kmer$kmer: $_\n" ;
			#}

		}
		close(VELG) ;

	} #end of determining kmer choice




		# skip if no file is produced
		unless (-s  "velvet_$kmerchoice.auto/$contig_L.$contig_R/contigs.fa") {
			print OUT "kmer choice:$kmerchoice: Zero contigs were produced!\n\n" ;	
			next ;
		}

		# check the file size and skip if file has size 0
		if ( -z "velvet_$kmerchoice.auto/$contig_L.$contig_R/contigs.fa") {
			print OUT "kmer choice:$kmerchoice: Zero contigs were produced!\n\n" ;
			next ;
		}



		# check the velvet scripts
		system("$checkvelvet velvet_$kmerchoice.auto/$contig_L.$contig_R/contigs.fa $kmerchoice velvet_$kmerchoice.auto/$contig_L.$contig_R/newcontigs.fa");
	
	
		# ssaha stage

		#ssaha the contigs
		#note: contigs are put into per file, because ssaha has the problem of not finding one or another..
		print "\nAligning the velvet contig...\n\n" ;

	
		#getting the contig in fasta format, for ssaha
		#note: contigs are put into per file, because ssaha has the problem of not finding one or another..
		if ( $contig_L ne 'xxxxxx' ) {
			open OUT_TMP, '>', "tmp/$PI/$PI.$random_no.contigs.L.fa" or die "cannot open tmp contig file" ;
			print OUT_TMP "$contigs_fasta{$contig_L}" ;
			close(OUT_TMP) ;
		}
	
		if ( $contig_R ne 'xxxxxx' ) {
			open OUT_TMP, '>', "tmp/$PI/$PI.$random_no.contigs.R.fa" or die "cannot open tmp contig file" ;
			print OUT_TMP "$contigs_fasta{$contig_R}" ;
			close(OUT_TMP) ;
		}


		my @output = () ;
		my @output_R = () ;

		if ( $contig_L ne 'xxxxxx' ) {
			@output = split /\n/ , `$ssaha -output cigar -best 1 tmp/$PI/$PI.$random_no.contigs.L.fa velvet_$kmerchoice.auto/$contig_L.$contig_R/newcontigs.fa | grep 'cigar'` ;
		}
		if ( $contig_R ne 'xxxxxx' ) {
			@output_R = split /\n/ , `$ssaha -output cigar -best 1 tmp/$PI/$PI.$random_no.contigs.R.fa velvet_$kmerchoice.auto/$contig_L.$contig_R/newcontigs.fa | grep 'cigar'` ;
		}

		push (@output, @output_R) ;
		#print "$kmer: @output" ;
		if (@output == 0) {

		}

		foreach (@output) {
			push( @{$velvet_contig_result{"$kmerchoice"}}, $_ )  ;
		}


	

	if ( $maxcontig_len < 250 ) {
		print OUT "misc: No contigs > 250bp was found!\n\n" ;
		next ;
	}


	#print out result
	my @final_ssaha_output = () ;
	if ( $velvet_contig_result{"$kmerchoice"} ) {
	    @final_ssaha_output = @{ $velvet_contig_result{"$kmerchoice"} } ;
	}
	else {
	    print OUT "misc: contig produced but no ssaha2 match to the original assembly!\n\n" ;
	    next ;
	}


			my $mapped_but_outside = 1 ;

			foreach ( @final_ssaha_output  ) {
				my @cigar_split = split /\s+/;
		
				# if the contig is within specified size; 5' end
				if ( $cigar_split[5] eq $contig_R && $cigar_split[6] <= $contig_ins_start{$contig_R}  ) {
					print OUT "$cigar_split[0] " ;
					print OUT "$contig_L.$contig_R.$cigar_split[1] " ;
					print OUT "@cigar_split[2..$#cigar_split]\n" ;

					unless ( $Illumina_contigs{"$contig_L.$contig_R.$cigar_split[1]"} ) {
						my $contig_fasta_seq = `awk '/$cigar_split[1]/ {getline; print }' velvet_$kmerchoice.auto/$contig_L.$contig_R/newcontigs.fa ` ;
						print OUTFA ">$contig_L.$contig_R.$cigar_split[1]\n" ;
						print OUTFA "$contig_fasta_seq" ;			
						$Illumina_contigs{"$contig_L.$contig_R.$cigar_split[1]"}++ ;

						$mapped_but_outside = 0 ;

					}

				}
				# if the contig is within specified size; 3' end
				elsif ( $cigar_split[5] eq $contig_L && $cigar_split[7] >= $contig_ins_end{$contig_L} ) {
					print OUT "$cigar_split[0] " ;
					print OUT "$contig_L.$contig_R.$cigar_split[1] " ;
					print OUT "@cigar_split[2..$#cigar_split]\n" ;

					unless ( $Illumina_contigs{"$contig_L.$contig_R.$cigar_split[1]"} ) {
						my $contig_fasta_seq = `awk '/$cigar_split[1]/ {getline; print }' velvet_$kmerchoice.auto/$contig_L.$contig_R/newcontigs.fa ` ;
						print OUTFA ">$contig_L.$contig_R.$cigar_split[1]\n" ;
						print OUTFA "$contig_fasta_seq" ;			
						$Illumina_contigs{"$contig_L.$contig_R.$cigar_split[1]"}++ ;

						$mapped_but_outside = 0 ;
					}

				}
				else {
					#print OUT "$_\toutside_range\n" ;
				}
	
			}


	print OUT "misc: contig produced but all mapped within contig ends!\n" if $mapped_but_outside == 1;
	print OUT "\n" ;

}


#last;
#$tmp_SC++ ;
#last if $tmp_SC == 3 ;

} #end of foreach loop

close(OUT) ;
close(OUTFA) ;


