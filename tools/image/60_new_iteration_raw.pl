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
if ($0 =~ /(\S+)\/60_new_iteration_raw.pl/ ) {
	#print "The path is $1\n" ; 	
	$path = $1 ; 
}


if (@ARGV != 10 ) {
	print "60_new_iteration.pl \n\n" ;
	print "has to be at the top directory\n\n" ;
	exit ;
}

my $dir_prefix = shift ;
my $new = shift ;
my $old = shift ;
my $prefix = shift ;
my $iteration = shift ;
my $iteration_all = shift ;
my $kmer = shift ;

#my $expcov = shift ;
#my $covcut = shift ;

my $toignore = shift ;
my $overhang = shift ;
my $velvet_insert_size = shift ; 

if ( $iteration_all != 0 ) {

mkdir "$new" or die "can not create new iteration dir $new !\n" ;
chdir "$new" ;


system("ln -s ../$old/new.read.placed image.read.placed") ;
system("ln -s ../$old/new.insert_size image.insert_size") ;
system("ln -s ../$old/new.fa image.contigs.fa") ;



print "all done! starting the new iteration at dir: $new\n" ;

print "commands to execute:\n $path/image.pl -dir_prefix $dir_prefix -prefix $prefix -iteration $iteration -all_iteration $iteration_all -toignore $toignore -overhang $overhang -kmer $kmer -vel_ins_len $velvet_insert_size\n\n" ;
system("$path/image.pl -dir_prefix $dir_prefix -prefix $prefix -iteration $iteration -all_iteration $iteration_all -toignore $toignore -overhang $overhang -kmer $kmer -vel_ins_len $velvet_insert_size") ;



}
else {
    print "\n\nall iteration have finished!\n\n" ;

}


