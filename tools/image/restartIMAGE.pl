#!/usr/bin/perl
#

#	author: jit
#
#
#
#

use strict;
use Getopt::Long;


my $path = '' ; 
if ($0 =~ /(\S+)\/restartIMAGE.pl/ ) {
	#print "The path is $1\n" ; 	
	$path = $1 ; 
}



if (@ARGV != 4 ) {
	print "\n\nIMAGE: restart image with different kmer\n" ;
	print "Author: jit\@sanger.ac.uk\n\n" ;

	print "\tUsage: $0 last_iteration_folder kmer number_of_iterations Illumina_lane_prefix\n" ;
	print "\t\tHas to be at the top directory\n" ;

	print "\n\tExample: $0 Ite_10 41 5 partitioned\n" ; 
	print "\t\tThis will run IMAGE by creating new iteration folder starting from Ite_11, for another five iterations of kmer 41, using partitioned_1.fastq and partitioned_2.fastq \n" ; 
	
	print "\tAdvanced usage: $0 -smalt_minscore 30 Ite_10 41 5 partitioned\n" ; 
	
	exit ;
}

my $SMALT_minScore = 0 ;

GetOptions (	                'smalt_minscore=i' => \$SMALT_minScore, 
) or die "Incorrect usage!\n";



my $dir_prefix = shift ;
my $kmer = shift ;
my $iteration_all = shift ;
my $illumina_prefix = shift ; 





if ( -d "$dir_prefix") {

}
else {
    print "folder $dir_prefix not present... Exiting...\n" ; 
    exit ; 
}




my $new = '' ;
my $old =  $dir_prefix;
my $prefix = '' ;
my $iteration = '' ;

if ( $dir_prefix =~ /(\D+)(\d+)/ ) {
	$prefix = $1 ; 
	$new = "$1" . ($2 + 1) ; 
	$iteration = $2 + 1 ;
}


my $toignore = 10 ;
my $overhang = 50 ;


if ( $iteration_all != 0 ) {

mkdir "$new" or die "can not create new iteration dir $new !\n" ;
chdir "$new" ;


system("ln -s ../$old/new.read.placed image.read.placed") ;
system("ln -s ../$old/new.insert_size image.insert_size") ;
system("ln -s ../$old/new.fa image.contigs.fa") ;



print "all done! starting the new iteration at dir: $new\n" ;

print "commands to execute:\n $path/image.pl -dir_prefix $prefix -prefix $illumina_prefix -iteration $iteration -all_iteration $iteration_all -toignore $toignore -overhang $overhang -kmer $kmer -smalt_minscore $SMALT_minScore\n\n" ;
system("$path/image.pl -dir_prefix $prefix -prefix $illumina_prefix -iteration $iteration -all_iteration $iteration_all -toignore $toignore -overhang $overhang -kmer $kmer -smalt_minscore $SMALT_minScore") ;



}
else {
    print "\n\nall iteration have finished!\n\n" ;

}


