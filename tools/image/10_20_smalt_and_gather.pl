#! /usr/bin/perl -w
#
#
# Copyright (C) 2009 by Pathogene Group, Sanger Center
#
# Author: JIT
# Description: 
#		a script that map the Solexa reads back to the reference/contigs
#		and parition them based on either ends of the contigs
#
#
#



use strict;
use warnings;



if ( @ARGV != 6) {
	print "$0 mate Flane Rlane insertsize output min_score\n\n" ;
	exit;
}

my $path = '' ;
if ($0 =~ /(\S+)\/10_20_smalt_and_gather.pl/ ) {
        #print "The path is $1\n" ;                                                                                                                                                
        $path = $1 ;
}


my $mate = shift ;
my $Flane = shift ;
my $Rlane = shift ;
my $range2 = shift ;
my $output = shift ;


# set minimum score
my $min_score = shift ;



my $range1 = 1   ;
my $insert_size = $range2 ;
my $contig_len_file = "../contigs.fa.len.txt" ;
my $insert = "../image.insert_size" ;


open READ_F, ">", "../partitioned_1.fastq" or die "oooops\n" ; 
open READ_R , ">", "../partitioned_2.fastq" or die "oooops\n" ; 



###############################################################################################
#read contig length file

my %contig_ins_start = ();
my %contig_ins_end = ();
my @contigs ;
my %contigs_len ;

my %contig_half_length = ();
my %contig_full_length = ();
open( IN, "$contig_len_file" ) or die "Cannot open $contig_len_file\n";

while (<IN>) {
    if (/(^\S+)\s+(\d+)/) {
	my $contig = $1 ;
	my $length = $2 ;
	push(@contigs, $contig) ;

	$contigs_len{$contig} = $length ;
	$contig_full_length{$contig} = $length ;
	$contig_half_length{$contig} = $length * 0.5 ;


	#$contig_ins_start{$contig} = $insert_size ;
	#$contig_ins_end{$contig} = $contigs_len{$contig} - $insert_size ;
    }

}
close(IN);


#print "using the file instead!\n" ;
	open( IN, "$insert" ) or die "Cannot open insert file\n";
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

###############################################################################################



#run ssaha....
#my $smalt_command = "/nfs/users/nfs_j/jit/bin/smalt-0.5.2/smalt map -x -i $range2 -f samsoft -o $output ref.fa $Flane $Rlane" ;
#my $smalt_command = " ~mh12/bin/smalt-0.5.4/smalt_x86_64 map -x -i $range2 -f samsoft -o $output ref.fa $Flane $Rlane" ;
#my $smalt_command = "/nfs/users/nfs_h/hp3/sw/arch/x86_64-linux/bin/smalt-0.5.5  map -x -i $range2 -f samsoft -o $output ref.fa $Flane $Rlane" ;


my $smalt_command = "$path/smalt_x86_64 map -x -i $range2 -f samsoft -o $output ref.fa $Flane $Rlane" ;


if ( $min_score != 0 ) {
    $smalt_command = "$path/smalt_x86_64 map -m $min_score -x -i $range2 -f samsoft -o $output ref.fa $Flane $Rlane" ;
}


print "executing smalt...\n" ;
system("$smalt_command");


print "parsing output...." ;
open( IN, "$output" ) or die "Cannot open $output (sam output)\n";

my $batch = 0 ;
if ($output =~ /(\d+)\.sam/) {
	print "batch output found\n" ;
	$batch = $1 ;
}


my %contig_reads ;

#my $count = 0 ;
my $total_reads = 0 ;

while (<IN>) {

    next if /^\@/ ; 


	my @first_read = split /\s+/, $_ ;
	my @second_read = split /\s+/ , <IN> ;

	my $F = $first_read[1] ;
	my $R = $second_read[1] ;

	my $partition = 0 ; 

	#print "@first_read\n" ;
	#print "@second_read\n" ;

	# flags of all for first pair...

	# if both unmapped
	if ($F == 77) {
	    $partition++ ; 
	}
	# if second read is unmapped and the first one is mapped
	if ($F == 73 || $F == 89) {
		if ( $first_read[3] >= $contig_ins_end{$first_read[2]}  ) {
			$contig_reads{$first_read[2]}{'Rend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;
		}
		elsif ( $first_read[3] <= $contig_ins_start{$first_read[2]} ) {
			$contig_reads{$first_read[2]}{'Lend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;
		}
		


		$partition++ ; 

	}
	# if first read is unmapped and the second read is mapped
	elsif ($R == 153 || $R == 137) {
		if ( $second_read[3] >= $contig_ins_end{$second_read[2]}  ) {
			$contig_reads{$first_read[2]}{'Rend'}{$second_read[0]} = "$first_read[9]\t$second_read[9]" ;
		}
		elsif ( $second_read[3] <= $contig_ins_start{$second_read[2]} ) {
			$contig_reads{$second_read[2]}{'Lend'}{$second_read[0]} = "$first_read[9]\t$second_read[9]" ;
		}

		$partition++ ; 
	}
	# mapped uniquely and within insert size
	elsif ($F == 83 || $F == 99) {
		#print "length: $contig_ins_end{$first_read[2]}\n" ;
		if ( $first_read[3] >= $contig_ins_end{$first_read[2]} && $first_read[7] >= $contig_ins_end{$first_read[2]}  ) {
			$contig_reads{$first_read[2]}{'Rend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;

			$partition++ ; 
		}
		elsif ( $first_read[3] <= $contig_ins_start{$first_read[2]} && $first_read[7] <= $contig_ins_start{$first_read[2]} ) {
			$contig_reads{$first_read[2]}{'Lend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;

			$partition++ ; 
		}
	}
	# mapped uniquely but outside specified insert size...
	elsif ($F == 81 || $F == 97) {
		if ( $first_read[3] >= $contig_ins_end{$first_read[2]} && $first_read[7] >= $contig_ins_end{$first_read[2]}  ) {
			$contig_reads{$first_read[2]}{'Rend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;

			$partition++ ; 
		}
		elsif ( $first_read[3] <= $contig_ins_start{$first_read[2]} && $first_read[7] <= $contig_ins_start{$first_read[2]} ) {
			$contig_reads{$first_read[2]}{'Lend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;

			$partition++ ; 
		}
	}
	# mapped uniquely but wrong orientation
	elsif ($F == 112 || $F == 199) {
		# not included

	    $partition++ ; 
	}
	# mapped uniquely in different contigs
	elsif ($F ==  65 || $F == 113) {
		# not included

		if ( $first_read[3] >= $contig_ins_end{$first_read[2]} || $first_read[3] <= $contig_ins_start{$first_read[2]}  ) {
			if ( $second_read[3] >= $contig_ins_end{$second_read[2]} || $second_read[3] <= $contig_ins_start{$first_read[2]}  ) {

				# just put in the first one for now....

				my $random = int(rand(2)) ;
			
				#print "$random\n" ;

				if ($random == 0 ) {
					if ( $first_read[3] >= $contig_ins_end{$first_read[2]} ) {
						$contig_reads{$first_read[2]}{'Rend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;
					}
					else {
						$contig_reads{$first_read[2]}{'Lend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;
					}
				}
				else {
					if ( $second_read[3] >= $contig_ins_end{$second_read[2]} ) {
						$contig_reads{$second_read[2]}{'Rend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;
					}
					else {
						$contig_reads{$second_read[2]}{'Lend'}{$first_read[0]} = "$first_read[9]\t$second_read[9]" ;
					}
				}

			

				#print "@first_read\t$contigs_len{$first_read[2]}\n" ;
				#print "@second_read\t$contigs_len{$second_read[2]}\n" ;	
				#$count++ ;		
				#last if $count == 100;
			}

		}


		$partition++ ; 
	}


	if ( $partition == 1 ) {
	    print READ_F '@' . "$first_read[0]\n$first_read[9]\n+\n$first_read[10]\n" ;
	    print READ_R '@' . "$second_read[0]\n$second_read[9]\n+\n$second_read[10]\n" ;
	}
	


}


print "done sorting..\n" ;

my @ends = qw / Lend Rend / ;
my $reads_used = 0 ;
my $out_dir = "../bridge/" ;
#mkdir $out_dir or die "can not make $out_dir\n" ;


#print out contigs for each file

open OUT, ">" , "$out_dir/$batch.list" or die $! ;

if ( $mate ne '0') {
    open OUT_MATE, ">", "$out_dir/$batch.matepairlist" or die $! ;
}


if ( $mate eq '0' ) {

    foreach my $contig (@contigs) {
	
	foreach my $end (@ends) {
	    
	    if ( $contig_reads{$contig}{$end} ) {			
		
		for my $readname (keys %{ $contig_reads{$contig}{$end} } ) {
		    my @reads = split /\t/ , $contig_reads{$contig}{$end}{$readname} ;
		    print OUT "$contig\t$end\t$readname\t$reads[0]\t$reads[1]\n" ;
		    $reads_used++ ;
		}			
		}
	}
    }
    
}
else {
    
    foreach my $contig (@contigs) {
	
        foreach my $end (@ends) {

            if ( $contig_reads{$contig}{$end} ) {
		
                for my $readname (keys %{ $contig_reads{$contig}{$end} } ) {
                    my @reads = split /\t/ , $contig_reads{$contig}{$end}{$readname} ;

		    if ( $readname =~ /$mate/ ) {
			print OUT_MATE "$contig\t$end\t$readname\t$reads[0]\t$reads[1]\n" ;
		    }
		    else {
			print OUT "$contig\t$end\t$readname\t$reads[0]\t$reads[1]\n" ;
		    }

                    $reads_used++ ;
                }
	    }
        }
    }

}


print "$reads_used pairs of reads are used to assemble gap regions\n" ;




