package MLocarna::SparseMatrix;

############################################################
##
## package of functions for handling 2D and 4D sparse matrices that
## are implemented as hashs of hashs
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
read_sparsematrix_2D
read_sparsematrix_4D
write_sparsematrix_2D
write_sparsematrix_4D
transpose_sparsematrix_2D
transpose_sparsematrix_4D
divide_sparsematrix_2D
divide_sparsematrix_2D_inplace
divide_sparsematrix_4D
scale_sparsematrix_4D
add_sparsematrix_2D
add_sparsematrix_2D_inplace
add_sparsematrix_4D
filter_sparsematrix_2D
filter_sparsematrix_4D
);

our %EXPORT_TAGS = ();

# your exported package globals go here,
# as well as any optionally exported functions
our @EXPORT_OK   = qw(
);


############################################################
## operations on 2- and 4-dimensional sparse matrices that are
## represented as hashs of hashs
## 

########################################
## add two 2-dimensional sparse matrices
sub add_sparsematrix_2D($$) {
    my ($m1,$m2) = @_;
        
    my %m3 = %{ $m1 };

    foreach my $i (keys %$m2) {
	foreach my $j (keys %{ $m2->{$i} } ) {
	    $m3{$i}{$j} += $m2->{$i}->{$j};
	}
    }
    
    return %m3;
}

########################################
## add_sparsematrix_2D_inplace($m1,$m2)
##
## Add a sparse matrices in-place
##
## @post $m1 = $m1 + $m2
##
########################################
sub add_sparsematrix_2D_inplace($$) {
    my ($m1,$m2) = @_;

    foreach my $i (keys %$m2) {
	foreach my $j (keys %{ $m2->{$i} } ) {
	    $m1->{$i}{$j} += $m2->{$i}->{$j};
	}
    }
}



########################################
## add two 4-dimensional sparse matrices
sub add_sparsematrix_4D($$) {
    my ($m1_ref,$m2_ref) = @_;

    my %m1 = %{ $m1_ref };
    my %m2 = %{ $m2_ref };

    my %m3 = %m1;


    foreach my $i (keys %m2) {
	foreach my $j (keys %{ $m2{$i} }) {
	    foreach my $k (keys %{ $m2{$i}{$j} }) {
		foreach my $l (keys %{ $m2{$i}{$j}{$k} }) {
		    $m3{$i}{$j}{$k}{$l} += $m2{$i}{$j}{$k}{$l};
		}
	    }
	}
    }
    
    return %m3;
}

sub divide_sparsematrix_2D($$) {
    my ($m_ref,$divisor) = @_;

    my %m = %{ $m_ref };
    my %r;

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    $r{$i}{$j} = $m{$i}{$j} / $divisor; 
	}
    }
    
    return %r;    
}

########################################
## divide_sparsematrix_2D_inplace($m,$divisor)
##
## @post $m /= $divisor
##
########################################
sub divide_sparsematrix_2D_inplace($$) {
    my ($m,$divisor) = @_;

    foreach my $i (keys %$m) {
	foreach my $j (keys %{ $m->{$i} }) {
	    $m->{$i}{$j} = $m->{$i}{$j} / $divisor; 
	}
    }
}


sub divide_sparsematrix_4D($$) {
    my ($m_ref,$divisor) = @_;

    my %m = %{ $m_ref };
    my %r;


    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    foreach my $k (keys %{ $m{$i}{$j} }) {
		foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
		    $r{$i}{$j}{$k}{$l} = $m{$i}{$j}{$k}{$l} / $divisor; 
		}
	    }
	}
    }

    return %r;    
}


########################################
## filter_sparsematrix_2D($m,$threshold)
##
## filter 2-dim matrix by probability threshold
##
## return result sparsematrix
########################################
sub filter_sparsematrix_2D($$) {
    my ($m_ref,$threshold) = @_;

    my %m = %{ $m_ref };
    my %r;

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    if ( $m{$i}{$j} >= $threshold ) {
		$r{$i}{$j} = $m{$i}{$j};
	    }
	}
    }
    
    return %r;    
}

########################################
## filter_sparsematrix_2D($m,$threshold)
##
## filter 4-dim matrix by probability threshold
##
## return result sparsematrix
########################################
sub filter_sparsematrix_4D($$) {
    my ($m_ref,$threshold) = @_;

    my %m = %{ $m_ref };
    my %r;

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    foreach my $k (keys %{ $m{$i}{$j} }) {
		foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
		    if ( $m{$i}{$j}{$k}{$l} >= $threshold ) {
			$r{$i}{$j}{$k}{$l} = $m{$i}{$j}{$k}{$l};
		    }
		}
	    }
	}
    }

    return %r;
}


sub scale_sparsematrix_2D($$) {
    my ($m_ref,$scale) = @_;

    my %m = %{ $m_ref };
    my %r;

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    $r{$i}{$j} = $m{$i}{$j} * $scale; 
	}
    }
    
    return %r;    
}

sub scale_sparsematrix_4D($$) {
    my ($m_ref,$scale) = @_;

    my %m = %{ $m_ref };
    my %r;


    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    foreach my $k (keys %{ $m{$i}{$j} }) {
		foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
		    $r{$i}{$j}{$k}{$l} = $m{$i}{$j}{$k}{$l} * $scale; 
		}
	    }
	}
    }

    return %r;    
}


# read from a sparse matrix file as written by locarna --write-match-probs 
sub read_sparsematrix_2D($) {
    my ($file) = @_;

    local *SM_IN;
    
    open(SM_IN,$file);

    my %h;
    
    while( my $line=<SM_IN> ) {
	if ( $line =~ /(\d+) (\d+) ([\d.e+-]+)/ ) {
	    $h{$1}{$2} = $3;
	}
    }
    
    close SM_IN;

    return %h;
}

# read from a sparse matrix file as written by locarna --write-match-probs 
sub read_sparsematrix_4D($) {
    my ($file) = @_;

    local *SM_IN;
    
    open(SM_IN,$file);

    my %h;
    
    while( my $line=<SM_IN> ) {
	if ( $line =~ /(\d+) (\d+) (\d+) (\d+) ([\d.e+-]+)/ ) {
	    $h{$1}{$2}{$3}{$4} = $5;
	}
    }
    
    close SM_IN;

    return %h;
}


sub write_sparsematrix_2D($$) {
    my ($m_ref,$file)=@_;
    my %m = %{ $m_ref };

    local *SM_OUT;
    
    open(SM_OUT,">$file") || die "Cannot write to $file."; 

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    print SM_OUT "$i $j $m{$i}{$j}\n"; 
	}
    }
    
    close SM_OUT;
}

sub print_sparsematrix_2D($$) {
    my ($m_ref)=@_;
    my %m = %{ $m_ref };

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    print "$i $j $m{$i}{$j}\n"; 
	}
    }
}


sub write_sparsematrix_4D($$) {
    my ($m_ref,$file)=@_;
    my %m = %{ $m_ref };

    local *SM_OUT;
    
    open(SM_OUT,">$file") || die "Cannot write to $file."; 

    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    foreach my $k (keys %{ $m{$i}{$j} }) {
		foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
		    print SM_OUT "$i $j $k $l $m{$i}{$j}{$k}{$l}\n"; 
		}
	    }
	}
    }
    
    close SM_OUT;
}


# transpose a sparse matrix implemented as hash of hashs
sub transpose_sparsematrix_2D($) {
    my ($m_ref) = @_;

    my %m = %{ $m_ref };
    
    my %r;
    
    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    $r{$j}{$i} = $m{$i}{$j}; 
	}
    }
    
    return %r;
}

# transpose a sparse matrix implemented as hash of hashs
sub transpose_sparsematrix_4D($) {
    my ($m_ref) = @_;

    my %m = %{ $m_ref };
    
    my %r;
    
    foreach my $i (keys %m) {
	foreach my $j (keys %{ $m{$i} }) {
	    foreach my $k (keys %{ $m{$i}{$j} }) {
		foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
		    $r{$k}{$l}{$i}{$j} = $m{$i}{$j}{$k}{$l};
		}
	    }
	}
    }
    
    return %r;
}

## ------------------------------------------------------------
1;
