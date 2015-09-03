#!/usr/bin/perl
#
#	to incorporate the information from velvet contig (in .afg)
#		to original ace file (which has been converted in .afg)
#
#	author: jit
#
#
#
#

use strict;

my $path = '' ; 
if ($0 =~ /(\S+)\/image_run_summary.pl/ ) {
	#print "The path is $1\n" ; 	
	$path = $1 ; 
}


if (@ARGV != 1) {
	print "\n\nIMAGE: Collect summaries from iteration directories\n" ;
	print "Author: jit\@sanger.ac.uk\n\n" ;

	print "\tUsage: $0 iteration_prefix \n" ; 
	print "\t\tHas to be at the top directory\n" ;

	print "\n\tExample: $0 Ite \n" ; 
	print "\t\tThis will collect statistics from Ite_1, Ite_2, Ite_3, ... etc\n\n\n" ; 
	
	exit ;
}

my $dir_prefix = shift ;



my @files = <*> ; 


print "The prefix is : $dir_prefix\n" ; 
print "iteration\tStarting_gaps\tGap_closed\tGap_extend_oneside\tGap_extend_bothside\n" ; 

my %result = () ; 

foreach my $file (@files) {

#    print "$file\n" ; 

    if ( $file =~ /^$dir_prefix(\d+)$/ ) {
#	print "here!\n" ; 

	my $iteration = $1 ; 

	my $totalgaps = `grep 'Total gaps' $file/walk2.summary | awk '{print \$5}'` ; chomp($totalgaps) ;
	my $gapclosed = `grep 'gap closed' $file/walk2.summary | awk '{print \$3}'` ; chomp($gapclosed) ; 
	my $gaponeside = `grep 'one side' $file/walk2.summary | awk '{print \$3}'` ;  chomp($gaponeside) ; 
	my $gapbothside = `grep 'both side' $file/walk2.summary | awk '{print \$3}'` ; chomp($gapbothside); 
	    
	#print "$file\t$totalgaps\t$gapclosed\t$gaponeside\t$gapbothside\n" ; 

	$result{$iteration} = "$totalgaps\t$gapclosed\t$gaponeside\t$gapbothside" ;
    }

}


for my $ite (sort  { $a <=> $b } keys %result ) {
    print "$ite\t$result{$ite}\n" ; 
}


