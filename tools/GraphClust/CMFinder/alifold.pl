#!/usr/bin/perl -w
use warnings;
use strict;
use POSIX;
use Getopt::Long;

my $usage = << "JUS";
  USAGE:  ./mloc2stockhol.pl -file STRING
  OPTIONS:
    -file		[MANDATORY]	File to convert into stockholm format


JUS

my $file;
my $opt_h = 0;
usage()

  unless GetOptions( "file=s" => \$file );

if ($opt_h) {
    print STDERR $usage;
    exit;
}

#######################
### CHECKLIST
#######################
if ( !$file ) {
    print "No file to convert....\n";
    print STDERR $usage;
    exit;
}

my $aln_file = $file;
my $tmp_path = "var/tmp";

## alifold result to get consensus structure string for infernal and some nice pictures
my $tmp_dir = "$tmp_path/alifold_$$";
my $currDir = getcwd;
mkdir($tmp_dir);
chdir($tmp_dir);
my @call_alifold =
  readpipe("RNAalifold -r --noLP --color --aln $aln_file 2>/dev/null");
my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
chomp($aln_cons_struct);
$aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
$aln_cons_struct = $1;
open( CONS, ">$aln_file.alifold" );
print CONS $call_alifold[0];
print CONS "$aln_cons_struct\n";
print CONS $2;
close(CONS);
system("mv alirna.ps $aln_file.alirna.ps");
system("mv aln.ps $aln_file.ps");
chdir($currDir);
