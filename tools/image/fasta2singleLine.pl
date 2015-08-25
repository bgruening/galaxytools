#!/usr/bin/perl

# Martin Hunt 02/10/08
# mh12@sanger.ac.uk

# Takes a fasta file and changes it so each sequence is on one line

use strict;
use warnings;

if ($#ARGV < 1){
    print "usage: fasta2singleLine.pl <input fasta file> <output fasta file> [1]\n",
          "where including the 1 will replace newlines with a space\n",
          "(useful for running on a qual file)\n\n";
    exit;
}

my $infile = $ARGV[0];
my $outfile = $ARGV[1];

# put whole file into an array, 1 sequence per element
open F, $infile or die ("Error opening $infile");
my $whole_file = join ("", <F>);    
close F;
chomp $whole_file;
$whole_file = substr ($whole_file,1);  # removes the > at the start of the file
my @sequences = split("\n>", $whole_file);


# write the output file
open F, ">$outfile" or die ("Error opening $outfile");

foreach (@sequences) {
    my @seq = split "\n";        # seq[0]=header, seq[1],seq[2],...=sequence
    print F ">$seq[0]\n";        # print seq header to file
    shift @seq;                  # take the header out of the array


    # put the sequence into 1 string (no newlines)
    my $out;



    if ($ARGV[2]){
        $out = join " ", @seq;
    }
    else{
        $out = join "", @seq;     
    }

    $out =~ s/ //g;

    print F "$out\n";            # print sequence to file
}

close F;
