package MLocarna::MatchProbs;

############################################################
##
## package for handling match probabilties and
## reliabilities. Comprises reading and writing, transformations,
## scoring of alignments, ...  Also contains closely related functions
## like max_weight_structure and write_dotplot
##
##
############################################################


use 5.008003;
use strict;
use warnings;

use MLocarna;
use MLocarna::Aux;
use MLocarna::SparseMatrix;

require Exporter;
    
# set the version for version checking
our $VERSION     = 1.00;

our @ISA         = qw(Exporter);
our @EXPORT      = qw(
read_bm_probs
write_bm_probs
read_am_probs
write_am_probs
write_bm_reliabilities
write_am_reliabilities
write_reliability_bars

average_basematch_probs
average_arcmatch_probs

consistency_transform_bm
consistency_transform_am

compute_amreliabilities_single_seq
compute_bmreliabilities_single_seq
compute_reliability

max_weight_structure
write_dotplot
write_dotplot_extended

evaluate_alignment
aln_reliability

);

our %EXPORT_TAGS = ();

# your exported package globals go here,
# as well as any optionally exported functions
our @EXPORT_OK   = qw(
$min_bm_prob
$min_am_prob
$scale_after_ct
);

our $min_bm_prob = 0.0005; # minimal prob of base match (used for pf alignment output)
our $min_am_prob = 0.0005; # minimal prob of arc match (used for pf alignment output)

our $scale_after_ct=0; # whether to scale probabilities after consistency transformation,
                      # scaling compensates for ignoring transition via gaps in the trafo


########################################
## sum_paired_prob_bm($bmprobs)
##
## $bmprobs reference to a hash of base match probabilities for all sequence pairs
##
## returns sum of all probabilities in the hash
##
########################################
sub sum_paired_prob_bm {
    my ($bmprobs) = @_;
    
    
    my @name_pairs = keys %$bmprobs;
    
    my $sum=0;
    
    foreach my $name_pair (@name_pairs) {
	foreach my $i ( keys %{ $bmprobs->{$name_pair} } ) {
	    foreach my $j ( keys %{ $bmprobs->{$name_pair}{$i} } ) {
		$sum += $bmprobs->{$name_pair}{$i}{$j};
	    }
	}
    }
    
    return $sum;
}

 ########################################
## sum_paired_prob_am
##
## $amprobs ref of hash of arc match probabilities
##
## returns sum of all probabilities in the hash
##
########################################
sub sum_paired_prob_am {
    my ($amprobs) = @_;
        
    my @name_pairs = keys %$amprobs;
    
    my $sum=0;
    
    foreach my $name_pair (@name_pairs) {
	foreach my $i ( keys %{ $amprobs->{$name_pair} } ) {
	    foreach my $j ( keys %{ $amprobs->{$name_pair}{$i} } ) {
		foreach my $k ( keys %{ $amprobs->{$name_pair}{$i}{$j} } ) {
		    foreach my $l ( keys %{ $amprobs->{$name_pair}{$i}{$j}{$k} } ) {
			$sum += $amprobs->{$name_pair}{$i}{$j}{$k}{$l};
		    }
		}
	    }
	}
    }

    return $sum;
}

########################################
## average_basematch_probs($alnA ref of alignment hash,$alnB ref of alignment hash,
##                         $bmprobs ref of base match probabilities hash)
##
## calculate the average base match probabilities for
## the alignment of the alignments alnA and alnB
##
## returns hash of averaged base match probabilities 
##
########################################
sub average_basematch_probs($$$) {
    my ($alnA,$alnB,$bmprobs) = @_;
    
    my @namesA = aln_names($alnA);
    my @namesB = aln_names($alnB);
    
    my %r; # result matrix

    for my $nameA (@namesA) {
	for my $nameB (@namesB) {
	    # projecterei

	    my @posmapA = project_seq($alnA->{$nameA});
	    my @posmapB = project_seq($alnB->{$nameB});
	    
	    my  @bmprobkeys=keys %$bmprobs;

	    my %m = %{ $bmprobs->{MLocarna::nnamepair($nameA,$nameB)} };

	    # compute sum over all sequence pairs
	    
	    foreach my $i (keys %m) {
		foreach my $j (keys %{ $m{$i} }) {
		    $r{$posmapA[$i]}{$posmapB[$j]} += $m{$i}{$j};
		}
	    }
	}
    }

    my $nm = ($#namesA+1)*($#namesB+1);
    
    return divide_sparsematrix_2D(\%r,$nm);
}


########################################
## average_arcmatch_probs($alnA ref of alignment hash,$alnB ref of alignment hash,
##                        $amprobs ref of arc match probabilities hash)
##
## calculate the average arc match probabilities for
## the alignment of the alignments alnA and alnB
##
## returns hash of averaged arc match probabilities 
##
########################################
sub average_arcmatch_probs($$$) {
    my ($alnA,$alnB,$amprobs) = @_;

    my @namesA = aln_names($alnA);
    my @namesB = aln_names($alnB);

    my %r; # result matrix

    for my $nameA (@namesA) {
	for my $nameB (@namesB) {
	    # projecterei
	    
	    my @posmapA = project_seq($alnA->{$nameA});
	    my @posmapB = project_seq($alnB->{$nameB});
	    
	    my %m = %{ $amprobs->{MLocarna::nnamepair($nameA,$nameB)} };
	    
	    # averagerei
	    
	    foreach my $i (keys %m) {
		foreach my $j (keys %{ $m{$i} }) {
		    foreach my $k (keys %{ $m{$i}{$j} }) {
			foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
			    $r{$posmapA[$i]}{$posmapA[$j]}{$posmapB[$k]}{$posmapB[$l]} += $m{$i}{$j}{$k}{$l};
			}
		    }
		}
	    }
	}
    }
	    
    my $nm = ($#namesA+1)*($#namesB+1);
    
    return divide_sparsematrix_4D(\%r,$nm);
}



########################################
## consistency_transform_single_am($mCA ref of hash,$mCB ref of hash)
##
## Performs consistency transformation on arc match probabilities for
## a single triple of sequences A, B, and C.  It re-estimates the
## probabilities for matches between A and B by match probabilities to
## a third sequence C.
##
## $mCA hash of arc match probabilities between C and A
## $mCB hash of arc match probabilities between C and B
##
## returns ref of hash %mAB of re-estimated match probabilities between A and B
##
########################################
sub consistency_transform_single_am($$) {
    my ($mCA,$mCB) = @_;

    my %mAB;
        
    foreach my $p (keys %$mCA) {
	foreach my $q (keys %{ $mCA->{$p} } ) {
	
	    if (exists $mCB->{$p} && exists $mCB->{$p}{$q}) {
		
		my %hA = %{ $mCA->{$p}{$q} };
		my %hB = %{ $mCB->{$p}{$q} };
		
		foreach my $i (keys %hA) {
		    foreach my $j (keys %{ $hA{$i} }) {
			foreach my $k (keys %hB) {
			    foreach my $l (keys %{ $hB{$k} }) {
				$mAB{$i}{$j}{$k}{$l} += $hA{$i}{$j} * $hB{$k}{$l}
			    }
			}
		    }
		}
	    }
	}
    }
    return \%mAB;
}

########################################
## consistency_transform_single_bm($mCA ref of hash,$mCB ref of hash)
##
## Performs consistency transformation on base match probabilities for
## a single triple of sequences A, B, and C.  It re-estimates the
## probabilities for matches between A and B by match probabilities to
## a third sequence C.
##
## $mCA hash of arc match probabilities between C and A
## $mCB hash of arc match probabilities between C and B
##
## returns ref of hash %mAB of re-estimated match probabilities between A and B
##
########################################
sub consistency_transform_single_bm($$) {
    my ($mCA_ref,$mCB_ref) = @_;

    my %mAB;
    
    my %mCA = %{ $mCA_ref };
    my %mCB = %{ $mCB_ref };
    
    # ATTENTION: actually we need to traverse the cut set of lists keys mCA and keys mCB
    foreach my $k (keys %mCA) {
	if (exists $mCB{$k}) {
	    my %hA = %{ $mCA{$k} };
	    my %hB = %{ $mCB{$k} };
	    
	    foreach my $i (keys %hA) {
		foreach my $j (keys %hB) {
		    $mAB{$i}{$j} += $hA{$i} * $hB{$j}
		}
	    }
	}
    }
    return \%mAB;
}

########################################
## consistency_transform_bm( $h_ref, $names_ref)
##
## get a hash indexed by name pairs of "sparse matrices" = match probs
## for all sequence pairs and the list of sequence names @names
## returns the transformed hash
##
## $h_ref reference to hash of base match probabilities
## $names_ref reference to list of names
##
########################################
sub consistency_transform_bm {
    my ( $h, $names ) = @_;
    
    my $num_seqs=@$names;

    my $sum_paired_prob;

    if ($scale_after_ct) {
       $sum_paired_prob = sum_paired_prob_bm($h);
    }
    
    # print_k_dim_hash(\%h,4,"");
    
    my %th; ## transformed hash
    
    my @name_pairs = keys %$h;
    
    foreach my $name_pair (@name_pairs) {
	my ($nameA,$nameB) = dhp($name_pair);
	if (exists $th{chp($nameB,$nameA)}) {
	    my %subhash = transpose_sparsematrix_2D( $th{chp($nameB,$nameA)} );
	    $th{$name_pair} = \%subhash;
	} else {
	    my %matrixAB;
	    
	    foreach my $nameC (@$names) {
		
		# careful, names in the hash are normalized!
		my $nnameC=MLocarna::get_normalized_seqname($nameC);
		
		# transform via $nameC and add to matrix AB
		
		my $matrixAB_add;
		
		if ($nnameC eq $nameA || $nnameC eq $nameB) {
		    $matrixAB_add = $h->{$name_pair}
		} else {
		    $matrixAB_add = 
			consistency_transform_single_bm( $h->{chp($nnameC,$nameA)},
							 $h->{chp($nnameC,$nameB)} );
		}
		
		add_sparsematrix_2D_inplace ( \%matrixAB, $matrixAB_add );
	    }
	
	    %matrixAB = filter_sparsematrix_2D( \%matrixAB, $num_seqs * $min_bm_prob);
	    
	    divide_sparsematrix_2D_inplace( \%matrixAB, $num_seqs );
	    
	    $th{$name_pair} = \%matrixAB;
	}
    }
    
    if ($scale_after_ct) {
	my $paired_prob_ratio = $sum_paired_prob/sum_paired_prob_bm(\%th);
	# printmsg 1, "Sum paired prob bm ratio: $paired_prob_ratio\n";
	
	foreach my $name_pair (@name_pairs) {
	    my %subhash = scale_sparsematrix_2D($th{$name_pair},$paired_prob_ratio);
	    $th{$name_pair} = \%subhash;
	}
    }
    
    return %th;
}

########################################
## read_bm_probs($dir directory name)
##
## read basematch probabilities from all files
## in directory of the form nameA-nameB
##
## $dir directory with base match probability files
##
## returns ref of hash of base match probabilities
##
########################################
sub read_bm_probs($) {
    my ($dir)=@_;
    
    my %bmprobs;
    
    local *MYIN;

    local *DIR;
    opendir(DIR, $dir);
    my @files = readdir DIR;
    closedir DIR;
        
    foreach my $file (@files) {
	if ($file =~ /(.+)-(.+)/) {
	    my $nameA=$1;
	    my $nameB=$2;
	    open(MYIN,"$dir/$file") || die "Cannot read $dir/$file";
	    while (my $line = <MYIN>) {
		if ($line =~ /(\d+) (\d+) ([\d.e+-]+)/) {
		    my $i=$1;
		    my $j=$2;
		    my $p=$3;
		    $bmprobs{MLocarna::chp($nameA,$nameB)}{$i}{$j}=$p;
		}
	    }
	    close MYIN;
	}
    }
    
    return \%bmprobs;
}


## write the basematch probabilities to files for later inspection by a user
sub write_bm_probs($$) {
    my ($dir,$bmprobs_ref)=@_;
    
    my %bmprobs = %{ $bmprobs_ref };


    local *MYOUT;
    
    mkdir "$dir" || die "Cannot make dir $dir";
    
    my @name_pairs = keys %bmprobs;

    for my $name_pair (@name_pairs) {
	my ($nameA,$nameB) = MLocarna::dhp($name_pair);

	my $outfile="$dir/$nameA-$nameB";
	open(MYOUT,">$outfile") || die "Cannot write to $outfile";
	
	my %h = %{ $bmprobs{$name_pair} };
	
	for my $i (keys %h) {
	    my %h2 = %{ $h{$i} };
	    for my $j (keys %h2) {
		printf MYOUT "$i $j %.8f\n", $h2{$j};
	    }
	}
	
	close MYOUT;
    }
}

########################################
## read_am_probs($dir directory name)
##
## read arc match probabilities from all files
## in directory of the form nameA-nameB
##
## $dir directory with arc match probability files
##
## returns ref of hash of arc match probabilities
##
########################################
sub read_am_probs($) {
    my ($dir)=@_;
  
    my %amprobs;
    
    local *MYIN;

    local *DIR;
    opendir(DIR, $dir);
    my @files = readdir DIR;
    closedir DIR;
    
    
    foreach my $file (@files) {
	if ($file =~ /(.*)-(.*)/) {
	    my $nameA=$1;
	    my $nameB=$2;
	    open(MYIN,"$dir/$file") || die "Cannot read $dir/$file";
	    while (my $line = <MYIN>) {
		if ($line =~ /(\d+) (\d+) (\d+) (\d+) ([\d.e+-]+)/) {
		    $amprobs{MLocarna::chp($nameA,$nameB)}{$1}{$2}{$3}{$4}=$5;
		}
	    }
	    close MYIN;
	}
    }
    
    return \%amprobs;
}


## write the basematch probabilities to files for later inspection by a user
sub write_am_probs($$) {
    my ($dir,$amprobs_ref)=@_;
    
    my %amprobs = %{ $amprobs_ref };
    
    local *MYOUT;
    
    mkdir "$dir" || die "Cannot make dir $dir"; 

    my @name_pairs = keys %amprobs;

    for my $name_pair (@name_pairs) {
	my ($nameA,$nameB) = MLocarna::dhp($name_pair);
	
	my $outfile="$dir/$nameA-$nameB";
	open(MYOUT,">$outfile") || die "Cannot write to $outfile";
	
	my %m = %{ $amprobs{$name_pair} };
		
	foreach my $i (keys %m) {
	    foreach my $j (keys %{ $m{$i} }) {
		foreach my $k (keys %{ $m{$i}{$j} }) {
		    foreach my $l (keys %{ $m{$i}{$j}{$k} }) {
			printf MYOUT "$i $j $k $l %.8f\n",$m{$i}{$j}{$k}{$l};
		    }
		}
	    }
	}
	close MYOUT;
    }
}

########################################
## compute_bmreliabilities_single_seq($name,\%aln,\%bmprobs,\%amprobs)
##
## compute bmreliabilities for a single sequence
##
## @returns pair \@res_bmrel_seq,\@res_bmrel_str
##
########################################
sub compute_bmreliabilities_single_seq {
    my ($nameA,$aln,$bmprobs,$amprobs) = @_;
    
    my @names = aln_names($aln);
    
    if ($#names < 0) { die "compute bmreliability single seq: empty alignment."; }

    my $aln_length=length($aln->{$names[0]});
        
    my @res_bmrel_seq; # result vector, bm contributions to bm rel
    ## accumulate sum of pairs of bmprobs in column $i in $res_bmrel_seq[$i]
    
    my @res_bmrel_str; # result vector, am contributions to bm rel
    ## accumulate sum of pairs of amprobs in column $i to some $j in $res_bmrel_str[$i]
    
    my %col2pos = %{ alicol2seqpos(\@names,$aln) }; ## mapping of alignment column to sequence position
    
    my %posmaps;
    
    for my $name (@names) {
	my @posmap=project_seq($aln->{$name});
	$posmaps{$name} = [ @posmap ];
    }

    ## initialize all $res_bmrel_seq[$i] and $res_bmrel_str[$i]
    for (my $i=1; $i<=$aln_length; $i++) {
	$res_bmrel_seq[$i] = 0;
	$res_bmrel_str[$i] = 0;
    }

    for (my $i=1; $i<=$aln_length; $i++) {
	
	if (substr($aln->{$nameA},$i-1,1) !~ /[A-Za-z]/ ) { 
	    $res_bmrel_seq[$i]  = 0;
	    $res_bmrel_str[$i] = 0;
	    next;
	}
	
	my @posmapA = @{ $posmaps{$nameA} };
	    
	for my $nameB (@names) {
	    if ($nameA eq $nameB) { next; }
		
	    if (substr($aln->{$nameB},$i-1,1) !~ /[A-Za-z]/ ) { next; }
	    
	    my $pA=$col2pos{$nameA}[$i];
	    my $pB=$col2pos{$nameB}[$i];
	    
	    ## 1.) contribution from base match
	    if (exists $bmprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pB}) {
		$res_bmrel_seq[$i] += $bmprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pB};
	    }
	    
	    ## 2.) contribution from arc matchs
	    
	    my @posmapB = @{ $posmaps{$nameB} };
	    
	    for my $pA2 ( keys %{ $amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA} } ) {
		
		my $i2 = $posmapA[$pA2];
		
		if ($i2<$i) { next; }
		
		if ( substr($aln->{$nameA},$i2-1,1) !~ /[A-Za-z]/ ) { next; }
		
		if ( substr($aln->{$nameB},$i2-1,1) !~ /[A-Za-z]/ ) { next; }
		
		my $pB2 =  $col2pos{$nameB}[$i2];
		
		if ( ! exists( $amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pA2}{$pB}{$pB2} ) ) { next; }
		
		my $add_prob=$amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pA2}{$pB}{$pB2};
		
		$res_bmrel_str[$i]  += $add_prob;
		$res_bmrel_str[$i2] += $add_prob;
	    }
	}
    }
    
    my $numpairs = $#names;
    
    for (my $i=1; $i<=$aln_length; $i++) {
	$res_bmrel_seq[$i] /= $numpairs;
	$res_bmrel_str[$i] /= $numpairs;
    }
    
    return (\@res_bmrel_seq,\@res_bmrel_str);
}

########################################
## compute_amreliabilities_single_seq($name sequence name, $aln ref of alingnment hash,
##                                $pairbrobs ref to pairprobs hash)
##
## compute arc match reliabilities for the given sequence
## in the given multiple alignment
##
## @param $name sequence name
## @param $aln ref of the mulitple alignment
## @param $amprobs ref of the hash of arc match probabilities of seqs in the alignment
## @param $pairprobs ref to hash of pair probabilities of the sequence
##
## NOTE: positions in the returned profile are alignment columns (not
## positions of the unaligned sequence)
##   ! alignment column indices are 1..'alignment length' as sequence positions are
##
## return ref of hash with arc match reliabilities
## 
########################################
sub compute_amreliabilities_single_seq($$$$) {
    
    my ($nameA,$aln,$amprobs,$pairprobs) = @_;
    
    my @names = aln_names($aln);

    my $aln_length=length($aln->{$names[0]});

    my %res_amrel; # result vector, am rel

    my %col2pos = %{ alicol2seqpos(\@names,$aln) }; ## mapping of alignment column to sequence position
    
    my @posmapA = project_seq($aln->{$nameA});
	    
    for (my $i=1; $i<=$aln_length; $i++) {
	## lookup in alignment string $aln{$nameA} at $i-1 (substr positions 0..strlen-1)
	if (substr($aln->{$nameA},$i-1,1) !~ /[A-Za-z]/ ) { next; } 
	
	for my $nameB (@names) {
	    if ($nameA eq $nameB) { next; }
		
	    if (substr($aln->{$nameB},$i-1,1) !~ /[A-Za-z]/ ) { next; }
		
	    my $pA=$col2pos{$nameA}[$i];
	    my $pB=$col2pos{$nameB}[$i];
	    
	    my @posmapB = project_seq($aln->{$nameB});
	    
	    for my $pA2 ( keys %{ $amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA} } ) {
		    
		my $i2 = $posmapA[$pA2];
		
		if ($i2<$i) { next; }
		
		if ( substr($aln->{$nameA},$i2,1) !~ /[A-Za-z]/ ) { next; }
		
		if ( substr($aln->{$nameB},$i2,1) !~ /[A-Za-z]/ ) { next; }
		
		my $pB2 =  $col2pos{$nameB}[$i2];
		
		if ( ! exists( $amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pA2}{$pB}{$pB2} ) ) { next; }
		
		my $add_prob=$amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pA2}{$pB}{$pB2};
		
		$res_amrel{"$i $i2"} += $add_prob; 
	    }
	}
    }
    
    ## now add the dot plot of the sequence alone
    
    foreach my $k (keys %$pairprobs) {
	$k =~ /(\S+) (\S+)/;
	my $p1=$1; my $p2=$2;
	
	my $i1=$posmapA[$p1];
	my $i2=$posmapA[$p2];
	
	$res_amrel{"$i1 $i2"} += $pairprobs->{$k};
    }
    
    ## finally average in res_amrel
    my $numpairs = @names ;  ## there are edges to k-1 many other
			     ## sequences and the pair probs from the
			     ## dot plot of the sequence
    
    for my $k (keys %res_amrel) {
	$res_amrel{$k} /= $numpairs;
    }
    
    return \%res_amrel;
}


########################################
## compute_reliability($aln ref to alignment hash,
##                     $bmprobs ref to base match probs hash,
##                     $amprobs ref to arc match probs hash)
##
## Compute reliabilities for the basematches in each column of a
## multiple alignment, as well as reliabilities for structural match
## of each two alignment columns.  Distinguish structure and sequence
## match contribution to base match reliabilities
##
## @returns triple \@res_bmrel_seq,\@res_bmrel_str, \%res_amrel
##
## alignment column indices are 1..'alignment length' for structures bmrel_seq,bmrel_str, amrel.
##
########################################
sub compute_reliability($$$) {
    my ($aln,$bmprobs,$amprobs) = @_;
    
    my @names = aln_names($aln);
    
    if ($#names < 0) { die "basematch_reliability: empty alignment."; }

    my $aln_length=length($aln->{$names[0]});
        
    my @res_bmrel_seq; # result vector, bm contributions to bm rel
    ## accumulate sum of pairs of bmprobs in column $i in $res_bmrel_seq[$i]
    
    my @res_bmrel_str; # result vector, am contributions to bm rel
    ## accumulate sum of pairs of amprobs in column $i to some $j in $res_bmrel_str[$i]
    
    my %res_amrel; # result vector, am rel
    
    my %col2pos = %{ alicol2seqpos(\@names,$aln) }; ## mapping of alignment column to sequence position
    
    my %posmaps;
    
    for my $name (@names) {
	my @posmap=project_seq($aln->{$name});
	$posmaps{$name} = [ @posmap ];
    }

    ## initialize all $res_bmrel_seq[$i] and $res_bmrel_str[$i]
    for (my $i=1; $i<=$aln_length; $i++) {
	$res_bmrel_seq[$i] = 0;
	$res_bmrel_str[$i] = 0;
    }

    for (my $i=1; $i<=$aln_length; $i++) {
	
	for my $nameA (@names) {
	    if (substr($aln->{$nameA},$i-1,1) !~ /[A-Za-z]/ ) { next; } 
	    my @posmapA = @{ $posmaps{$nameA} };
	    
	    for my $nameB (@names) {
		if ($nameA le $nameB) { next; }
		
		if (substr($aln->{$nameB},$i-1,1) !~ /[A-Za-z]/ ) { next; }
		
		my $pA=$col2pos{$nameA}[$i];
		my $pB=$col2pos{$nameB}[$i];
		
		## 1.) contribution from base match
		if (exists $bmprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pB}) {
		    $res_bmrel_seq[$i] += $bmprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pB};
		}
		
		## 2.) contribution from arc matchs
		
		my @posmapB = @{ $posmaps{$nameB} };
	    		
		for my $pA2 ( keys %{ $amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA} } ) {
		    
		    my $i2 = $posmapA[$pA2];
		    
		    if ($i2<$i) { next; }
		    
		    if ( substr($aln->{$nameA},$i2-1,1) !~ /[A-Za-z]/ ) { next; }
		    
		    if ( substr($aln->{$nameB},$i2-1,1) !~ /[A-Za-z]/ ) { next; }
		    
		    my $pB2 =  $col2pos{$nameB}[$i2];
		    
		    #printerr "MLocarna::nnamepair($nameA,$nameB) $pA $pA2 $pB $pB2\n";
		    if ( ! exists( $amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pA2}{$pB}{$pB2} ) ) { next; }
		    
		    my $add_prob=$amprobs->{MLocarna::nnamepair($nameA,$nameB)}{$pA}{$pA2}{$pB}{$pB2};
		    
		    $res_amrel{"$i $i2"} += $add_prob; 
		    
		    #print "add $nameA $nameB $i $i2 [ $pA $pA2 $pB $pB2 ] $add_prob\n";
		    $res_bmrel_str[$i]  += $add_prob;
		    $res_bmrel_str[$i2] += $add_prob;
		}
	    }
	}
    }
    
    
    my $numpairs = ($#names+1)*$#names/2;
    
    for (my $i=1; $i<=$aln_length; $i++) {
	$res_bmrel_seq[$i] /= $numpairs;
	$res_bmrel_str[$i] /= $numpairs;
    }
    
    for my $k (keys %res_amrel) {
	$res_amrel{$k} /= $numpairs;
    }

    return (\@res_bmrel_seq,\@res_bmrel_str, \%res_amrel);
}



########################################
## write_bm_reliabilities($file filename, $bmrels_seq ref of list, $bmrels_str ref of list)
##
## write reliabilities to files
##
## $file name of output file
## $bmrels_seq ref of list of sequence  match reliabilities per column
## $bmrels_str ref of list of structure match reliabilities per column
##
## alignment column indices are 1..'alignment length'
##
########################################
sub write_bm_reliabilities {    
    my ($file, $bmrels_seq, $bmrels_str) = @_;
    
    open(BMREL,">$file") || die "Cannot write to $file";
    my $aln_length=@$bmrels_seq-1; # (!)
    for (my $i=1; $i<=$aln_length; $i++) {
	print BMREL "$i $bmrels_seq->[$i] $bmrels_str->[$i]\n";
    }
    close BMREL;
}

########################################
## write arc match reliabilities to file
##
## implicitely, alignment column indices are 1..'alignment length'
## (same for intern representation and file)
##
########################################
sub write_am_reliabilities {    
    my ($file, $amrels) = @_;
    
    open(AMREL,">$file")  || die "Cannot write to $file";
    foreach my $k (keys %$amrels) {
	print AMREL "$k $amrels->{$k}\n";
    }
    close AMREL;
}

########################################
## write reliability plot to the screen
##
## note column indices 1..'alignment length'
##
sub write_reliability_bars {
    my ($bmrels_seq, $bmrels_str) = @_;
    
    ## compute a sum of pairs of the base match probabilities for each alignment column in result.aln
     
    my $reso=10;
    for (my $i=0;$i<$reso; $i++) {
	print sprintf("-%3d%%              ",100*($i+1)/$reso);
	my $aln_length=@$bmrels_seq-1; # (!)
	for (my $j=1; $j<=$aln_length; $j++) {
	    my $b=$bmrels_seq->[$j];
	    my $a=$bmrels_str->[$j];# * ($mea_beta/100);
	    
	    if ($a>$i/$reso) {
		print "#";
	    } elsif ($a+$b>$i/$reso) {
		print "*";
	    } else {
		print " ";
	    }
	}
	print "\n";
    }
}

########################################
## write_dotplot_extended($filename,$sequence,$pairprobinfo,$colors1,$colors2)
##
## write a postscript dotplot of am probs/reliabilities to a file
## postscript code is adapted from Vienna RNA package RNAfold -p output
##
## $filename write to this filename
## $sequence sequence string
## $pairprobinfo ref of hash of pair probability information
##   ! the keys of the hash are string "$i $j", the ends of the pair
##   ! indices in $pairprobs need to be in range 1..'sequence length'
##
## $colors1 color table 1 for upper right triangle
## $colors2 color table 2 for lower left triangle
########################################
sub write_dotplot_extended {
    my ($filename,$sequence,$pairprobs,$colors1,$colors2) = @_;
    
    my $data="";
    foreach my $k (keys %$pairprobs) {
	$k =~ /(\d+) (\d+)/ || die "write_dotplot_extended: wrong keys in pairprobs.";
	my $i=$1;
	my $j=$2;
	if ( $i < 1 ) {
	    printerr "write_dotplot_extended: Out of range i: $i\n";
	}
	if ( $j > length($sequence) ) {
	    printerr "write_dotplot_extended: Out of range j: $j\n";
	}

	
	my @pairprobinfo = split /\s+/,$pairprobs->{$k};
	
	my $p1=$pairprobinfo[0];
	my $p2=$pairprobinfo[1];
	my $w1=$pairprobinfo[2];
	my $w2=$pairprobinfo[3];

	my $colnum1=@$colors1;
	my $colnum2=@$colors2;
	
	if ($p1!=0) {
	    if (!defined($w1)) { # do we have a weight for p1
		$w1=0;
	    }
	    $data .= $colors1->[int($w1*($colnum1-1))]." setrgbcolor ";
	    $data .= "$i $j ".sqrt($p1)." ubox\n";
	}

	if (defined($p2) && ($p2!=0)) { # do we have a weight for p1
	    if (!defined($w2)) { # do we have a weight for p1
		$w2=0;
	    }
	    $data .= $colors2->[int($w2*($colnum2-1))]." setrgbcolor ";
	    $data .= "$i $j ".sqrt($p2)." lbox\n";
	}
    }
    my $pscode="%!PS-Adobe-3.0 EPSF-3.0
%%Title: RNA Dot Plot
%%Creator: MLocARNA
%%CreationDate: Wed May 13 11:16:08 2009
%%BoundingBox: 66 211 518 662
%%DocumentFonts: Helvetica
%%Pages: 1
%%EndComments

%Options: 
% 
%This file contains the square roots of the base pair probabilities in the form
% i  j  sqrt(p(i,j)) ubox

%%BeginProlog
/DPdict 100 dict def
DPdict begin
/logscale false def
/lpmin 1e-05 log def

/box { %size x y box - draws box centered on x,y
   2 index 0.5 mul sub            % x -= 0.5
   exch 2 index 0.5 mul sub exch  % y -= 0.5
   3 -1 roll dup rectfill
} bind def

/ubox {
   logscale {
      log dup add lpmin div 1 exch sub dup 0 lt { pop 0 } if
   } if
   3 1 roll
   exch len exch sub 1 add box
} bind def

/lbox {
   3 1 roll
   len exch sub 1 add box
} bind def

/drawseq {
% print sequence along all 4 sides
[ [0.7 -0.3 0 ]
  [0.7 0.7 len add 0]
  [-0.3 len sub -0.4 -90]
  [-0.3 len sub 0.7 len add -90]
] {
   gsave
    aload pop rotate translate
    0 1 len 1 sub {
     dup 0 moveto
     sequence exch 1 getinterval
     show
    } for
   grestore
  } forall
} bind def

/drawgrid{
  0.01 setlinewidth
  len log 0.9 sub cvi 10 exch exp  % grid spacing
  dup 1 gt {
     dup dup 20 div dup 2 array astore exch 40 div setdash
  } { [0.3 0.7] 0.1 setdash } ifelse
  0 exch len {
     dup dup
     0 moveto
     len lineto 
     dup
     len exch sub 0 exch moveto
     len exch len exch sub lineto
     stroke
  } for
  [] 0 setdash
  0.04 setlinewidth 
  currentdict /cutpoint known {
    cutpoint 1 sub
    dup dup -1 moveto len 1 add lineto
    len exch sub dup
    -1 exch moveto len 1 add exch lineto
    stroke
  } if
  0.5 neg dup translate
} bind def

end
%%EndProlog
DPdict begin
%delete next line to get rid of title
270 665 moveto /Helvetica findfont 14 scalefont setfont (dot.ps) show

/sequence { (\\
$sequence\\
) } def
/len { sequence length } bind def

72 216 translate
72 6 mul len 1 add div dup scale
/Helvetica findfont 0.95 scalefont setfont

drawseq
0.5 dup translate
% draw diagonal
0.04 setlinewidth
0 len moveto len 0 lineto stroke 

drawgrid

%data starts here
$data
showpage
end
%%EOF
";

    local *OUT;
    open(OUT,">$filename") || die "Cannot write to $filename";
    print OUT $pscode;
    close OUT;

    return;
}


########################################
## write_dotplot($filename,$sequence,$pairprobs)
##
## write a simple postscript dotplot of am probs/reliabilities to a file
## postscript code is adapted from Vienna RNA package RNAfold -p output
##
## $filename write to this filename
## $sequence sequence string
## $pairprobs ref of hash of pair probabilities
##   ! the keys of the hash are string "$i $j", the ends of the pair
##   ! indices in $pairprobs need to be in range 1..'sequence length'
##
########################################
sub write_dotplot($$$) {
    my ($filename,$sequence,$pairprobs) = @_;
    
    my @colors=("0 0 0");
    
    write_dotplot_extended($filename,$sequence,$pairprobs,\@colors,\@colors);
}


########################################
## max_weight_structure($len,$base_weights,$arc_weights,$struct_weight,$threshold)
##
## Computes the structure with maximal sum of weights
##
## $len length of sequence
##   ! all indices need to be in range 1..$len
## $base_weight ref of hash of base weights
## $arc_weight ref of hash of arc weights
## $struct_weight multiply arc weights by this factor 
##
## $threshold controls which base pairs appear as brackets in output
##
## implements a candidate list approach to speed up computation
##
## returns pair of score and structure (as dot-bracket string using
## different bracket symbols to distinguish base pairs of different 'strength')
##
########################################
sub max_weight_structure($$$$$) {
    my ($len,$base_weights,$arc_weights,$struct_weight,$threshold) = @_;
    
    my @Nmat; ## a Nussinov like DP-matrix
    my @Tmat; ## Traceback matrix

    my $d_min=4; ## minimal difference $j-$i for a basepair ($i,$j)
    
    
    ## Nussinov like initialization
    for (my $i=1; $i<=$len; $i++) {
	$Nmat[$i][$i-1]=0;
	for (my $k=$i; $k<=$len && $k<$i+$d_min; $k++) {
	    my $bw = (@$base_weights==0)?0:$base_weights->[$k];
	    
	    $Nmat[$i][$k] = $Nmat[$i][$k-1] + $bw;
	    $Tmat[$i][$k]=0;
	}
    }
    
    ## Do Nussinov algorithm with candidate list approach in order to use sparsity
    ##
    
    for (my $i=$len; $i>=1; $i--) {
	my @arc_ends=(); ## candidate list for right ends of arcs (i,k)
	
	for (my $j=$i+$d_min; $j<=$len; $j++) {
	    
	    my $bw = (@$base_weights==0)?0:$base_weights->[$i]; ## support 0 base weights by passing empty array of base weights
	    $Nmat[$i][$j] = $Nmat[$i+1][$j] + $bw;
	    $Tmat[$i][$j] = -1;
	    
	    foreach my $k (@arc_ends) {
		if (! defined $Nmat[$k+1][$j]) {
		    printerr "$i $k $j\n";
		}
		my $score = $Nmat[$i+1][$k-1] + $Nmat[$k+1][$j] + $arc_weights->{"$i $k"} * $struct_weight;
		
		if ($score > $Nmat[$i][$j]) {
		    $Nmat[$i][$j] = $score;
		    $Tmat[$i][$j] = $k;
		}
	    }
	    
	    if ((exists $arc_weights->{"$i $j"})) {
		my $score = $Nmat[$i+1][$j-1] + $arc_weights->{"$i $j"} * $struct_weight;
		
		if ($score > $Nmat[$i][$j]) {
		    $Nmat[$i][$j] = $score;
		    $Tmat[$i][$j] = $j;
		
		    push @arc_ends, $j; ## only if there is an arc (i,j) and the score is maximal for N[i][j], we need to consider j as a new candidate
		                        ## (this follows the "candidate list" idea of Mihal Ziv-Ukelson and co-authors)
		}
	    }
	}
    }



    # read total score from Nmat
    my $total_score=$Nmat[1][$len];
    
    # DEBUGGING:
    # print Nmat
    
#     print "\n";
#     for (my $a=1; $a<=$len; $a++) {
# 	for (my $b=1; $b<=$len; $b++) {
# 	    if ($a>=$b) {print "xxxx ";} else { 
# 		printf "%4.1f ",$Nmat[$a][$b];
# 	    }
# 	}
# 	print "\n";
#     }

    # compute trace back
    
    my @stack;
    push @stack, (1,$len);

    my @trace;
    $trace[0]="";
    for (my $i=1; $i<=$len; $i++) {
	$trace[$i]=".";
    }

    while ($#stack != -1) {
		    
	my $i=$stack[-2];
	my $j=$stack[-1];
	$#stack -= 2;
	
	if ($i+$d_min>$j) { next; }
	
	if ($Tmat[$i][$j]==-1) {
	    push @stack,($i+1,$j);
	} elsif ($Tmat[$i][$j]>0) {
	    my $k = $Tmat[$i][$j];
	    if  ($arc_weights->{"$i $k"}>=$threshold) {
		$trace[$i]="(";
		$trace[$k]=")";
	    } elsif ($arc_weights->{"$i $k"}>=$threshold/2) {
		$trace[$i]="{";
		$trace[$k]="}";
	    } else {
		$trace[$i]="`";
		$trace[$k]="'";
	    }

	    push @stack,($i+1,$k-1);
	    push @stack,($k+1,$j);
	}
    }
    
    return ($total_score,join('', @trace));
}

########################################
## consistency_transform_am( $h reference to hash, $names ref of list )
##
## get a hash indexed by name pairs of "sparse matrices" = arc match
## probs for all sequence pairs
##
## $h     ref of hash of arc match probabilities
## $names ref of list of names
##
## returns the transformed hash
########################################
sub consistency_transform_am {
    my ( $h, $names ) = @_;
    
    my $num_seqs=@$names;
    
    my $sum_paired_prob;
    if ($MLocarna::scale_after_ct) {
       $sum_paired_prob = sum_paired_prob_am($h);
    }

    my %th;
    
    my @name_pairs = keys %$h;
    
    foreach my $name_pair (@name_pairs) {
	my ($nameA,$nameB) = dhp($name_pair);
	
	if ($nameA eq $nameB) {
	    $th{$name_pair} = $h->{$name_pair};
	} else {
	    if (exists $th{chp($nameB,$nameA)}) {
		my %subhash = transpose_sparsematrix_4D( $th{chp($nameB,$nameA)} );
		$th{$name_pair} = \%subhash;
	    } else {
		
		my %matrixAB;
		
		foreach my $nameC (@$names) {

		    # careful, names in the hash are normalized!
		    my $nnameC=MLocarna::get_normalized_seqname($nameC);
		    
		    # transform via $nameC and add to matrix AB
		    
		    my $matrixAB_add;
		    
		    if ($nnameC eq $nameA || $nnameC eq $nameB) {
			$matrixAB_add = $h->{$name_pair};
		    } else {
			$matrixAB_add = consistency_transform_single_am( $h->{chp($nnameC,$nameA)} , $h->{chp($nnameC,$nameB)} );
		    }
		    
		    %matrixAB = add_sparsematrix_4D ( \%matrixAB, $matrixAB_add );
		}
		
		%matrixAB = filter_sparsematrix_4D( \%matrixAB, $num_seqs * $min_am_prob);
		
		%matrixAB = divide_sparsematrix_4D( \%matrixAB, $num_seqs );
		
		$th{$name_pair} = \%matrixAB;
	    }
	}
    }


    if ($MLocarna::scale_after_ct) {
	
	my $paired_prob_ratio = $sum_paired_prob/sum_paired_prob_am(\%th);
	#printmsg 1, "Sum paired prob am ratio: $paired_prob_ratio\n";
	
	foreach my $name_pair (@name_pairs) {
	    my %subhash = scale_sparsematrix_4D($th{$name_pair},$paired_prob_ratio);
	    $th{$name_pair} = \%subhash; 
	}
    }
    
    return %th;
}


############################################################
##
## internally used functions
##


########################################
## alicol2seqpos($names ref of names list, $aln ref of alignment)
##
## compute mapping from alignment column to sequence position for
## each sequence in a given multiple alignment
## @param $names_ref ref of list only compute for given names
## @param $aln_ref ref of hash, the multiple alignment
##
## @returns ref of mapping (hash of arrays)
##
## alignment cols are in [1..alilen]
## sequence positions are in [1..seqlen]
## the hash entries are only valid for non-gap positions
##
########################################
sub alicol2seqpos($$) {
    my ($names_ref, $aln_ref) = @_;
    my @names = @{ $names_ref };
    my $aln_length=length($aln_ref->{$names[0]});

    my %col2pos;

    for my $name (@names) {
	my @clist;
	my $c=0;
	for (my $i=0; $i<$aln_length; $i++) {
	    if ( substr($aln_ref->{$name},$i,1) =~ /[A-Za-z]/ ) {
		$c++;
	    }
	    $clist[$i+1]=$c;
	}
	$col2pos{$name}=\@clist;
    }
    return \%col2pos;
}


########################################
## evaluate_alignment($aln,$bmprobs,$amprobs,$output)
##
## Evaluate an alignment by its reliabilities according to given
## base match and arc match probabilities
##
## $aln     ref of  alignment as alignment hash 
## $bmprobs ref of base match probs hash
## $amprobs ref of arc match probs hash 
## $output  boolean flag, whether to print to screen
##
## returns pair of reliability score and most reliable structure
##
########################################
sub evaluate_alignment {
    my ($aln,$bmprobs,$amprobs,$output) = @_;

    my $len=aln_length($aln);
    
    my ($bmrels_seq, $bmrels_str, $amrels) 
	= compute_reliability($aln,$bmprobs,$amprobs);

    if ($output) {
	if (3 >= $MLocarna::Aux::verbosemode) {
	    write_reliability_bars($bmrels_seq, $bmrels_str);
	}
	printmsg 3, "\n";
    }
        
    my ($score,$rel_str)=max_weight_structure($len,$bmrels_seq,$amrels,3.0,0);
    
    $score /= $len;
    $score *= 100;
    
    
    return ($score,$rel_str);
}

########################################
## aln_reliability($aln,$bmprobs,$amprobs)
##
## compute reliability score of alignment $aln given $bmprobs and $amprobs
## as averaged column reliability (unweighted structure+sequence rel.) 
##
########################################
sub aln_reliability($$$) {
    my ($aln,$bmprobs,$amprobs) = @_;

    my ($bmrels_seq, $bmrels_str, $amrels) 
	= compute_reliability($aln,$bmprobs,$amprobs);

    
    my $score=0;
    
    my $len=aln_length($aln);
    
    for (my $i=1; $i<=$len; $i++) {
	my $b=$bmrels_seq->[$i];
	my $a=$bmrels_str->[$i];
	
	$score += $a+$b;
    }
    
    $score /= $len;
    $score *= 100;
    
    return $score;
}


## ------------------------------------------------------------
1;
