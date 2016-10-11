#! /usr/bin/perl -w

# Zizhen Yao
# CVS $Id: cmfinder.pl,v 3.1 2006/03/07 19:22:09 yzizhen Exp $

use Getopt::Long qw(:config no_ignore_case); 

########## TO DO ##############
#~ change the path to gloibal


#$path= ".";
#$path= "/home/sokhoyae/Desktop/GalaxyProject/galaxy/tools/CMFinder";
#$path= $ENV{CMfinder};
$blast_path=$ENV{BLAST};

#default parameters

$CAND=40;
$N=3;
$MAXSPAN=100;
$MINSPAN=30;
$FRACTION=0.8;
$MIN_HAIRPIN = 1;
$MAX_HAIRPIN = 1;
$NO_BLAST = 0;
$verbose = 0;
$help = 0;

if (!GetOptions(
	   "h" => \$help,
	   "v" => \$verbose,
	   "b" => \$NO_BLAST,
	   "c=i" => \$CAND,
	   "n=i" => \$N,
	   "m=i" => \$MINSPAN,
	   "M=i" => \$MAXSPAN,
	   "f=f" => \$FRACTION,
	   "s=i" => \$MIN_HAIRPIN,
	   "S=i" => \$MAX_HAIRPIN
	   )
    ){
    print STDERR "Invalid options\n";
    print_help();
    exit(1);
}

if ($help) {
    print_help();
    exit(0);
}

if ($MAX_HAIRPIN < $MIN_HAIRPIN){
    $MAX_HAIRPIN  = $MIN_HAIRPIN;
}

if ($MAX_HAIRPIN == $MIN_HAIRPIN) {
    $HAIRPIN=$MIN_HAIRPIN;
}
else{
    $HAIRPIN = "$MIN_HAIRPIN.$MAX_HAIRPIN";
}

if ($NO_BLAST == 0 && (! defined $blast_path || ! -e "$blast_path/blastn" || ! -e "$blast_path/xdformat")){
    print STDERR "Can not find BLAST. Search without BLAST\n";
    $NO_BLAST = 1;
}

if (scalar @ARGV==0){
    print STDERR "No sequence file is specfied\n";
    print_help();
    exit(1);
}

$SEQ= shift @ARGV;

if ($verbose) {
    print <<OPTION;
N=$N
SEQ=$SEQ 
CAND=$CAND 
MAXSPAN=$MAXSPAN 
MINSPAN=$MINSPAN 
FRACTION=$FRACTION 
MIN_HAIRPIN=$MIN_HAIRPIN 
MAX_HAIRPIN=$MAX_HAIRPIN 

OPTION
}



system("./candf -c $CAND -o $SEQ.h$HAIRPIN.cand -m $MINSPAN -M $MAXSPAN -s $MIN_HAIRPIN -S $MAX_HAIRPIN $SEQ");

if ($NO_BLAST == 0) {
    #build blast database
    system("$blast_path/xdformat -n $SEQ 2>/dev/null");
    system("$blast_path/blastn $SEQ $SEQ -notes -top -W 8 -noseqs> $SEQ.blast");
    system("$path/parse_blast.pl $SEQ.blast > $SEQ.match"); 		    
    system("$path/cands -n $N -f $FRACTION -m $SEQ.match $SEQ $SEQ.h$HAIRPIN.cand");
}
else{ 
	
    system("./cands -n $N -f $FRACTION $SEQ $SEQ.h$HAIRPIN.cand");
}

for($i=1; $i <= $N; $i++) {
    if (-s "$SEQ.h$HAIRPIN.cand.$i"){
	system("./canda $SEQ.h$HAIRPIN.cand.$i  $SEQ $SEQ.align.h$HAIRPIN.$i");
	system("./cmfinder -o $SEQ.motif.h$HAIRPIN.$i -a $SEQ.align.h$HAIRPIN.$i $SEQ $SEQ.cm.h$HAIRPIN.$i >> $SEQ.h$HAIRPIN.out.$i");
    }
}

if (! $verbose){
    system("rm $SEQ.*cand*");
    system("rm $SEQ.*align*");
    system("rm $SEQ.*out*");
    if ($NO_BLAST == 0) {
	system("rm $SEQ.blast");
	system("rm $SEQ.match");
	system("rm $SEQ.xn*");
    }
}

sub print_help{
    print STDERR <<HELP;
cmfinder.pl [options] SEQ
Options:
    -b               Do not use BLAST search to locate anchors
    -v               Verbose. Print running information, and save intermediate results.
    -c <number>      The maximum number of candidates in each sequence. Default 40. No bigger than 100.
    -m <number>      The minimum length of candidates. Default 30
    -M <number>      The maximum length of candidates. Default 100
    -n <number>      The maximum number of output motifs. Default 3
    -f <number>      The fraction of the sequences expected to contain the motif. Default 0.80
    -s <number>      The number of stem-loops in the motif
    -h               Show this list
HELP
}
