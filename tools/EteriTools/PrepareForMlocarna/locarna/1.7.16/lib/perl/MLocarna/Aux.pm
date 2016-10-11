package MLocarna::Aux;

############################################################
##
## package of auxilliary functions
##
## put functions here that are used by multiple packages and don't fit
## somewhere else.
##
############################################################

use 5.008003;
use strict;
use warnings;

require Exporter;
    
# set the version for version checking
our $VERSION     = 1.00;

our @ISA         = qw(Exporter);
our @EXPORT      = qw(
printmsg
printerr
systemverb
systemverb_withinput
dhp
chp
subtract_list
project_seq
is_gap
aln_length
aln_size
aln_names
consensus_sequence
);



our %EXPORT_TAGS = ();

# your exported package globals go here,
# as well as any optionally exported functions
our @EXPORT_OK   = qw(
$verbosemode
);


our $verbosemode=3;


########################################
## printmsg($verblevel level of verbosity, $message message string)
##
## print message depending on verbosity-level
## ## and write to protocol file
##
## ## level -1 = only protocol file
## level 0  = print to screen only moreverbose
## level 2  = print to screen for verbose
## level 3  = print to screen
sub printmsg($$) {
    my ($verblevel, $message) = @_;
    
    ## print LOG "$message";

    if ($verblevel>=$verbosemode) {
	print "$message";
    }
}

########################################
## print error message
sub printerr($) {
    my ($message) = @_;
    
    ## print LOG "$message";

    print STDERR "$message";
}


########################################
## make systemcall and print call if verbose-mode
sub systemverb($) {
    my ($cmd)=@_;
    printmsg 1,"$cmd\n";
    printmsg 1,readpipe("$cmd");
}

########################################
## make systemcall with input and print call if verbose-mode
sub systemverb_withinput($$) {
    my ($input,$cmd)=@_;
    printmsg 1,"+++$input+++ >>> $cmd\n";
    
    $cmd.=">/dev/null" unless $verbosemode>0;
    my $fh;
    open($fh,"|$cmd");
    print $fh $input;
    close $fh;
}


## compose hash pair
sub chp($$) {
    my ($nameA,$nameB)=@_;
    return "$nameA#$nameB";
}

## decompose hash pair
sub dhp($) {
    my ($name_pair)=@_;
    $name_pair =~ /([^#]*)#([^#]*)/;
    return ($1,$2);
}


########################################
## subtract_list($l1,$l2)
##
## return list containing all elements of $l1 that don't occur in $l2
##
sub subtract_list {
    my ($l1, $l2) = @_;

    my @res;
    
    foreach my $x (@$l1) {
	my $found=0;
	foreach my $y (@$l2) {
	    if ($x eq $y) {
		$found=1;
		last;
	    }
	}
	if ($found==0) { push @res, $x; }
    }
    return @res;
}

########################################
## is_gap($s)
## returns whether $s is a gap symbol (different to A-Za-z
##
sub is_gap($) {
    my ($s)=@_;
    return ($s !~ /^[A-Za-z]$/); 
}

########################################
## project_seq($alig_str)
##
## produce an array that yields for each sequence position the corresponding
## position in the given alignment string
## sequence positions in [1..seqlen], alignment positions in [1..alilen]
##
########################################
sub project_seq($) {
    my ($alig_str) = @_;
    
    my @posmap;
    
    my $len=length($alig_str);
    
    my $j=1;
    for (my $i=1; $i<=$len; $i++) {
	if ( ! is_gap(substr($alig_str,$i-1,1)) ) {
	    $posmap[$j]=$i;
	    $j++;
	}
    }
    
    return @posmap;
}

########################################
## return size of alignment, i.e. number of sequences in the alignment
## @param $aln ref to alignment hash
##
## ATTENTION: aln is not allowed to contain constraint extensions!
sub aln_size($) {
    my $aln = shift;
    
    my @ks = keys %$aln;

    #if (grep /#[S,C,LONG]$/,@ks) {
    #	print STDERR "WARNING: Potentially you discovered a bug in mlocarna. Please avoid sequence names that end in #S, #C or #LONG: @ks\n";
    #}
    
    return $#ks+1;
}

########################################
## return sequence names of the alignment
## @param $aln ref to alignment hash
##
## ATTENTION: aln is not allowed to contain constraint extensions!
sub aln_names($) {
    my $aln = shift;
    
    my @ks = keys %$aln;

    #if (grep /#[S,C,LONG]$/,@ks) {
    #	print STDERR "WARNING: Potentially you discovered a bug in mlocarna. Please avoid sequence names that end in #S, #C or #LONG: @ks\n";
    #}
    
    return @ks;
}

########################################
## return size of alignment, i.e. number of sequences in the alignment
## @param $aln ref to alignment hash
##
## ATTENTION: sequence names in aln are not allowed to contain symbols '#',
## since this symbol is reserved for the constraint tags
##
sub aln_size_with_constraints($) {
    my $aln = shift;
    
    my @ks = keys %$aln;
    @ks = grep !/#/,@ks;
    return $#ks+1;
}

########################################
## aln_length($aln ref of alignment)
##
## return length of alignment strings
## (assume that all strings have same length)
sub aln_length($) {
    my $aln = shift;
    my @ks = keys %$aln;
    return length( $aln->{$ks[0]} );
}


# return the consensus sequence of an alignment
# for each column take the most frequently occurring symbol
#
sub consensus_sequence {
    my ($aln_ref) = @_;
       
    my $len=aln_length($aln_ref);
    my %aln = %{ $aln_ref };
        
    # count occurence of symbols in each alignment column 
    # and determine consensus (for each column take best count)

    my $consensus="";
    
    for (my $col=0; $col<$len; $col++) {
	
	my %counts;
	
	foreach my $name (keys %aln) {
	    my $sym = substr $aln{$name},$col,1;
	    $counts{$sym}++;
	}
	
	my $best=-1;
	my $best_sym="_";
	foreach my $sym (keys %counts) {
	    if (($counts{$sym} > $best) || (($counts{$sym} == $best) && $sym eq "-")) {
		$best = $counts{$sym};
		$best_sym=$sym;
	    }
	}
	$consensus .= $best_sym;
    }
    
    return $consensus;
}


## ------------------------------------------------------------
1;
