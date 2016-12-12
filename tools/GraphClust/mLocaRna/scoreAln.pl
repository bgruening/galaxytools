#!/usr/bin/perl -w
# -*-Perl-*-
# Author: Dominic Rose
# scoreAln.pl

use FindBin;
use lib "$FindBin::Bin";
use Getopt::Long;
use strict;
use Data::Dumper;
use Math::Round;

#adapt to each program
my $usage = << "JUS";
  usage:    perl $0 -i ALIGNMENT_FILE -f FORMAT [CLUSTALW|TBA|MAF] -o [0|1]  
  
  options:  -i    An alignment file 	                  [REQUIRED]
            -f    The format of the alignment             [REQUIRED]
	          [CLUSTALW|TBA|MAF]
            -o    Detailed output (1) or no details, only the score (0)
	    -s	  Which score				  [REQUIRED]
	          [sop|mpi]
	    
  purpose:  Calculate sum-of-pairs score, mean pairwise identity for given alignments
  results:  For each alignment the score is written to STDOUT  
JUS

my ( $opt_i, $opt_f, $opt_o, $opt_s ) = ( "", "", 0, "" );

my $alignLength = 0;

# read sequences
my @filedata;
my $seq = '';
my $id  = '';
my $cnt = 0;
my @keys;
my $cl     = 0;
my $qtSeqs = 0;

usage()
  unless GetOptions(
    "i=s" => \$opt_i,
    "f=s" => \$opt_f,
    "o=i" => \$opt_o,
    "s=s" => \$opt_s
  );

if ( !$opt_f || !$opt_i ) {
    print STDERR $usage;
    exit;
}

# ------ Start ----------

my %seq;

main();

sub main {
    $alignLength = read_sequences();
    die "ERROR: Alignment length is invalid: '$alignLength' for aln='$opt_i'"
      if ( $alignLength <= 0 );
    $qtSeqs = keys %seq;

    # sort keys (only once)
    foreach my $k ( sort keys %seq ) {
        push @keys, $k;
    }

    # get columnLength once
    my $column = get_column( 0, @keys );
    $cl = length $column;

    if ( $opt_s =~ /sop/ ) {
        compute_sop();
    }
    elsif ( $opt_s =~ /mpi/ ) {
        compute_mpi();
    }
}    # end of main

# prints SUM_PF_PAIRS SCORE
sub compute_sop {
    my $qtOverallScore = 0;

    my $qtPairs = 0;    # = binomial coefficient
    $qtPairs = binc( $cl, 2 );

    # for each column
    for ( my $i = 0 ; $i < $alignLength ; $i++ ) {
        my $column = get_column( $i, @keys );
        print "\n$i\. column = '$column'\n" if ( $opt_o == 1 );
        my $score = calcSumOfPairsScore( $column, $cl, $qtPairs );
        print "score = '$score'\n" if ( $opt_o == 1 );
        $qtOverallScore += $score;
    }

    # normalize #seqs
    my $tmp = nearest( .01, ( 2 / ( $qtSeqs - 1 ) * $qtOverallScore ) );

    # normalize score to 0..1
    #my $tmp = $qtOverallScore / ($qtPairs*$alignLength);

    if ( $opt_o == 1 ) {
        print
"\n--  overall score = 2 / ($qtSeqs -1) * $qtOverallScore = $tmp  --\n\n";

       #print "\n--  overall score = $qtOverallScore / $qtPairs = $tmp  --\n\n";
    }
    else {
        print "$tmp\n";
    }
}

# prints MEAN PAIRWISE IDENTITY
sub compute_mpi {

    my @aln = ();

    foreach my $k (@keys) {

        #print $seq{$k}{originalSeq};
        push @aln, [ split( //, $seq{$k}{originalSeq} ) ];
    }

    #print Dumper(@aln);
    #exit;

    my $pairs   = 0;
    my $matches = 0;

    for my $i ( 0 .. $#aln ) {
        for my $j ( $i + 1 .. $#aln ) {
            for my $k ( 0 .. ( @{ $aln[0] } - 1 ) ) {
                if ( ( $aln[$i][$k] ne '-' ) or ( $aln[$j][$k] ne '-' ) ) {
                    if ( $aln[$i][$k] eq $aln[$j][$k] ) {
                        $matches++;
                    }
                    $pairs++;
                }
            }
        }
    }

    my $mpi = sprintf( "%.4f", $matches / $pairs );

    if ( $opt_o == 1 ) {
        print "MPI = $mpi\n";
    }
    else {
        print "$mpi\n";
    }

}

# ----------------------------------------------------------------
# calculate binomial coefficient
sub binc {
    my ( $n, $k ) = @_;
    return 0 if ( $k > $n );
    return 0 if ( $k < 0 );
    if ( $k > int( $n / 2 ) ) {
        $k = $n - $k;
    }

    my $rv = 1;
    for my $j ( 0 .. $k - 1 ) {
        $rv *= $n - $j;
        $rv /= $j + 1;
    }
    return $rv;
}

# ----------------------------------------------------------------
# See http://mathworld.wolfram.com/LucasCorrespondenceTheorem.html.
# Write n and k in base-p notation, with digits n_i and k_i.  Then (n choose k)
# is equivalent mod p to the product of the (n_i choose k_i)'s.

# produce all pairs of nucleotides for a given column of an alignment
# calculate and return sum-of-pairs score
sub calcSumOfPairsScore() {
    my ( $column, $cl, $qtPairs ) = @_;

    my $score = 0;

    for ( my $i = 0 ; $i < $cl ; $i++ ) {
        print "$i:\n" if ( $opt_o == 1 );

        for ( my $j = 0 ; $j < $i ; $j++ ) {
            print "$i-$j " if ( $opt_o == 1 );
            $score += getScoreOfPair( $column, $i, $j );
        }
        print "\n" if ( $opt_o == 1 );
    }

    #my $tmp = nearest(.01, $score/$qtPairs);
    #print "Normalizing: $score / $qtPairs = $tmp\n" if($qtPairs>0);

    #return $tmp;
    return $score;
}

# sore a pair
sub getScoreOfPair() {
    my ( $column, $i, $j ) = @_;

    my $s = 0;
    my ( $iC, $jC ) = ( '', '' );

    $iC = substr $column, $i, 1;
    $jC = substr $column, $j, 1;

    if ( $iC ne "-" && $iC ne "N" && $iC eq $jC ) {
        $s = 1;
    }

    print "\n$i=$iC eq $j=$jC -> s=$s\n" if ( $opt_o == 1 );

    return $s;
}

# PURPOSE: Gets the column of an alignment as a string
# PARAMS:  $c: Index of column [0,MAX]
sub get_column {
    my ( $c, @keys ) = @_;

    my $column = '';

    foreach my $k (@keys) {
        $column .= substr $seq{$k}{originalSeq}, $c, 1;
    }

    return $column;

}    # end of sub get_column

sub read_sequences {
    open( IN, "< $opt_i" )
      or die "ERROR (read sequences): Could not open file '$opt_i'\n\n";
    @filedata = <IN>;
    close(IN);

    %seq = ();

    if ( $opt_f =~ /CLUSTALW/ ) {
        readClustalw();
    }
    elsif ( $opt_f =~ /TBA/ ) {
        readTba();
    }
    elsif ( $opt_f =~ /MAF/ ) {
        readMAF();
    }
    else {
        die "ERROR: Wrong aln format!\n";
    }

    my @Keys = keys %seq;

    return length $seq{ $Keys[0] }{originalSeq};

}    # end of sub read_sequences

sub readClustalw() {
    foreach my $line (@filedata) {
        if ( $line =~ /^\s*CLUSTAL\sW\s*$/ ) {
            next;
        }
        elsif ( $line =~ /^\s*$/ ) {

            $cnt = 0
              ; # Alignments are spread over several lines. A blank line indicates
                # new line of multiple alignment!
            next;
        }
        elsif ( $line =~ /^\s*[\*\s]+\s*$/ ) {
            next;
        }
        elsif ( $line =~ /^\s*(\S*)\s+(\S*)\s*$/ ) {
            $cnt++;
            $id = $cnt . "."
              . $1;   # To make sure that IDs are unique as DIALIGN cuts off IDs

            if ( !exists $seq{$id} ) {
                $seq{$id}{originalSeq} = '';
                $seq{$id}{valid}       = 1;
            }

            my $subseq = $2;
            $subseq =~ s/\s//g;
            $seq{$id}{originalSeq} .= $subseq;
        }
    }
}

sub readTba() {
    for ( my $i = 0 ; $i <= $#filedata ; $i++ ) {
        my $line = $filedata[$i];

        if ( $line =~ /^#$/ ) {
            next;
        }

        # >Nsp.genome:605198:75:+:tba_score=9506.0
        elsif ( $line =~ /^\>(.+)$/ ) {
            $cnt++;
            my $id = $1;

            $i++;
            my $seq = $filedata[$i];
            chomp $seq;
            $seq{$id}{originalSeq} .= $seq;
        }
    }
}

sub readMAF() {
    for ( my $i = 0 ; $i <= $#filedata ; $i++ ) {
        my $line = $filedata[$i];

        if ( $line =~ /^#$/ ) {
            next;
        }

        #s       Eco_1   21079   102     +       1       GGTCCAACT
        elsif ( $line =~ /^s\t(\S+)\t(\d+)\t(\d+)\t([+|-])\t(\d+)\t(\S+)$/ ) {
            $cnt++;
            my $id = $1 . "_" . $2 . "_" . $3;

            $seq{$id}{originalSeq} = $6;
            chomp $seq{$id}{originalSeq};
        }
    }

    #print Dumper %seq;
}

