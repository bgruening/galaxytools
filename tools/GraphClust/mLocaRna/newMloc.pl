#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/ min max /;
use Cwd qw(abs_path getcwd);
our @EXPORT_OK = qw(
  $node_sym
);

our $node_sym = "\$\$nodesym";

my (
    $center_fa_file, $tree_file, $tree_matrix, $p,
    $max_diff_am,    $tau,       $max_diff,    $path,
    $data_map,       $plfold_minlen
) = @ARGV;

my $OPTS_locarna_paligs =
" -p $p --max-diff-am $max_diff_am --tau $tau --max-diff $max_diff --indel-open -400 --indel -200 --struct-weight 180 ";
my $OPTS_locarna_maligs =
" -p $p --max-diff-am $max_diff_am --tau $tau --max-diff $max_diff --alifold-consensus-dp ";
my $OPTS_locarna_p_model =
" -p $p --max-diff-am $max_diff_am --tau $tau --max-diff $max_diff --struct-weight 180 --plfold-span 150 --plfold-winsize 200 --temperature 180 --mea-beta 400 --consistency-transformation ";
my $OPTS_locarna_model =
" -p $p --max-diff-am $max_diff_am --tau $tau --max-diff $max_diff --alifold-consensus-dp --indel-open -400 --indel -200 --struct-weight 180 ";
my $center_subtree_min   = 3;
my $center_subtree_max   = 7;
my $center_tree_type     = 3;
my $center_model_type    = 5;
my $center_tree_aligs    = 1;
my $center_skip_unstable = 0;
my $vrna_path            = ".";
my $cmfinder_path        = ".";

my $tree_aligs_local = 0;
$tree_aligs_local = 1
  if ( $OPTS_locarna_paligs =~ /local/
    || $OPTS_locarna_paligs =~ /normalized/ );

my $dp_dir = "dp";

my $use_sim_stree = 0;
my $use_prob_alig = 0;
$use_prob_alig = 1 if ( $center_tree_aligs == 2 );

my $mloc_opts = $OPTS_locarna_maligs;
$mloc_opts = $OPTS_locarna_p_model if $use_prob_alig;

if ( !$use_prob_alig ) {

    my $rnafold_opts   = "--noLP ";
    my $rnaplfold_opts = " --noLP -c 0.0005 -L 200 -W 300";

    #  my $plfold_minlen  =  210 ;

    my $fold_opts =
"--vrna-path $vrna_path --fasta $center_fa_file --tgtdir $dp_dir --new-only ";
    $fold_opts .= "--switch-len $plfold_minlen ";
    $fold_opts .= "--rnafold-opts " . $rnafold_opts;
    $fold_opts .= "--rnaplfold-opts " . $rnaplfold_opts;

    #print "fold opts = $fold_opts\n";
    system("perl $path/foldFasta.pl $fold_opts");
}

computeTreeAligsLocP(
    "maligs.loc_new",    $center_fa_file,
    $tree_file,          $center_subtree_min,
    $center_subtree_max, $use_sim_stree,
    $use_prob_alig,      $mloc_opts,
    $dp_dir,             "SUBTREES"
);

my $in_eval_mode = 0;

my $cluster_id = "MISS";
my $tree_dir   = "Tree";
system("mkdir $tree_dir");

my $subtrees = rankSubtrees( "SUBTREES", "$tree_file", $in_eval_mode );
my $best_subtrees_file = "bestSubtrees";
open( BEST, ">$best_subtrees_file" );

my $rank           = 0;
my @subtrees_final = ();

foreach
  my $t ( sort { $b->{SCORE_SIMPLE} <=> $a->{SCORE_SIMPLE} } @{$subtrees} )
{

    $rank++;
    print BEST getNodeInfo( $t, $cluster_id, $rank ) . "\n";

    ## no model if subtree contains overlapping sequences (frags)
    next if ( $t->{OVERLAP} );

    ## skip unstructured subtrees
    next if ( $t->{MFE} >= 0 && $center_skip_unstable );

    #next if ( $t->{RNAZ}->[4] >= 0 && $center_skip_unstructured );

    ## otherwise potential model alignment
    push( @subtrees_final, $t );
}

close(BEST);

makeCol("$best_subtrees_file");

if ( !@subtrees_final ) {
    print
"\nNo subtree/cluster with model criteria found in tree! No model will be build. Exit...\n\n";
    exit;
}

print "Found "
  . scalar(@subtrees_final)
  . " models. Build model alignment now for best subtree... \n";
my $model_dir = "MODEL";
system("mkdir -p $model_dir");

my @fa_center = read_fasta_file($center_fa_file);
my $subset_fa = writeSubsetFasta( \@fa_center, $subtrees_final[0]->{NAMES},
    "$model_dir/model.tree.fa", 1 );

## default file name = $center_model_type == 1
my $tree_aln_file = $subtrees_final[0]->{FILE};
system("cp $tree_aln_file $model_dir/best_subtree.aln");
system("cp $tree_aln_file.ps $model_dir/best_subtree.aln.ps");

if ( $center_model_type == 2 ) {
    print "center_model_type = 2\n";
    ## realign with locarna and diff opts 'OPTS_locarna_model'
    print "\n...realign best subtree with locarna...\n";
    mlocarna_center( $center_fa_file, $model_dir, $dp_dir, 0 );
    $tree_aln_file = "$model_dir/results/result.aln";

}
elsif ($center_model_type == 3
    || $center_model_type == 4
    || $center_model_type == 5 )
{

    if ( $center_tree_aligs != 2 && !-e "$model_dir/results/result.aln" ) {
        ## realign with LocarnaP as subtrees were aligned with locarna (with $center_tree_aligs==1 we already have locarnaP alignments)
        print "\n...realign best subtree with locarnaP...\n";
        mlocarna_center( $subset_fa, $model_dir, $dp_dir, 1 );
        $tree_aln_file = "$model_dir/results/result.aln";
    }

}

## in case we have realigned best subtree we create the new alifold files (annotated colored postscript)
# if ( !( -e "$tree_aln_file.ps" && -e "$tree_aln_file.alifold" ) ) {
#     aln2alifold( $tree_aln_file, '.', $vrna_path );
#     print "mtnuma aln2alifold i mej \n";
# }

## default is best subtree alignment
#my $final_tree_prefix = "$model_dir/model.tree.aln";
# system("cp $tree_aln_file $final_tree_prefix");
# system("cp $tree_aln_file.ps $final_tree_prefix.ps");
# system("cp $tree_aln_file.alirna.ps $final_tree_prefix.alirna.ps");
# system("cp $tree_aln_file.alifold $final_tree_prefix.alifold");
#
# ## reliabillity signal
# if ( -e "$tree_aln_file.rel_signal" ) {
#   system("cp $tree_aln_file.rel_signal $final_tree_prefix.rel_signal");
#   system("cp $tree_aln_file.rel_plot.pdf $final_tree_prefix.rel_plot.pdf")
#    if ( -e "$tree_aln_file.rel_plot.pdf" ); ## in case R is not working correctly, plot is not there
# }

sub mlocarna_center {
    my $fasta    = $_[0];
    my $dir      = $_[1];
    my $dpDir    = $_[2];
    my $use_locP = $_[3];

    #  my $OPTS_locarna_p_model = $CONFIG{OPTS_locarna_p_model};
    #  my $OPTS_locarna_model   = $CONFIG{OPTS_locarna_model};
    #  my $loc_path             = $CONFIG{PATH_LOCARNA};

    #  my $vrna_path            = $CONFIG{PATH_VRNA};

    my $loc_pp_dir = "$dir/input";
    system("mkdir -p $loc_pp_dir");

    my @fa = read_fasta_file($fasta);
    foreach my $key ( keys %{ $fa[0] } ) {
        system("ln -f -s $dpDir/$key $loc_pp_dir/$key")
          if ( -e "$dpDir/$key" && !$use_locP );

    }


#system("mlocarna --sequ-local true --tgtdir $dir $fasta $in_mlocarna_opts");
#system( "mlocarna $OPTS_locarna_model --verbose --skip-pp --tgtdir $dir $fasta > $dir/locarna.out 2>&1 ");
#  my $currDir = getcwd;
#  print "curr dir before mlocaran = $currDir \n";
    system(
"mlocarna $OPTS_locarna_model  --skip-pp --verbose --tgtdir $dir $fasta > $dir/locarna.out 2>&1"
    );

    #}

    ## cleanup
#  system("rm -r -f $dir/input $dir/intermediates $dir/probs $dir/results/single_reliabilities $dir/results/ali* $dir/alifold.out");

}

sub read_fasta_file {
    my ( $file, $make_unique ) = @_;
    my $FUNCTION = "read_fasta_file in Sequences.pm";

    my $id         = "";
    my $seqstring  = "";
    my %fasta      = ();
    my %header     = ();
    my @order      = ();
    my $line       = "";
    my %meta       = ();
    my %seq_meta   = ();
    my $uniq_count = 0;

    open( IN_HANDLE, "<$file" )
      || die "ERROR in $FUNCTION: "
      . "Couldn't open the following file in package Tool,"
      . " sub read_fasta_file: $file\n";

    while ( $line = <IN_HANDLE> ) {
        chomp($line);

        # header (can contain one space after > symbol)
        if ( $line =~ /^\>\s?(\S+)\s*([\S*\s*]*)/ ) {
            if ($id) {
                if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {

                    #          $uniq_count++;
                    #          $id .= "_$uniq_count";
                    #          print "Warning! Make Seq-id unique! now $id\n";
                    die "ERROR in $FUNCTION: "
                      . "multiple sequence id '$id', consider using function "
                      . "read_fasta_with_nonunique_headers instead of read_fasta_file";
                }
                $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
                $fasta{$id} = $seqstring;
                $meta{$id}  = {%seq_meta};    ## anonymous hash reference

                $seqstring = "";
                undef(%seq_meta);
            }

            if ($make_unique) {
                $uniq_count++;
                $id = $1 . "_$uniq_count";
            }
            else {
                $id = $1;
            }

            my $head = $2;
            $head = "" if ( !$head );
            $header{$id} = $head;
            push( @order, $id );
        }
        elsif ( $line =~ /(.+)\s+(#\S+)\s*$/ && $id ) {

            if ( exists $seq_meta{$2} ) {
                $seq_meta{$2} .= $1;
            }
            else {
                $seq_meta{$2} = $1;
            }

        }
        else {
            $seqstring .= $line if ($id);
        }
    }

    if ($id) {
        if ( defined $fasta{$id} and ( $fasta{$id} ne $seqstring ) ) {

            #$uniq_count++;
            #$id .= "_$uniq_count";
            #print "Warning! Make Seq-id unique! now $id\n";
            die "ERROR in $FUNCTION: "
              . "multiple sequence id '$id', consider using function "
              . "read_fasta_with_nonunique_headers instead of read_fasta_file";
        }
        $seqstring =~ s/\s*//g;    ## do not allow spaces in sequence
        $fasta{$id} = $seqstring;
        $meta{$id}  = +{%seq_meta};    ## anonymous hash reference
        $seqstring  = "";
        undef(%seq_meta);

    }

    return ( \%fasta, \@order, \%header, \%meta );
}

## min_seqs : further evaluate subtrees of sim-tree with >$min_seqs leafs
## max_seqs: further evaluate subtrees of sim-tree with <$max_seqs leafs
## use_sim_tree: uses tree structure from sim tree or not for realigning subtree
sub computeTreeAligsLocP {
    my $tgtDir            = $_[0];
    my $fa_tree_file      = $_[1];
    my $tree_file         = $_[2];
    my $min_seqs          = $_[3];
    my $max_seqs          = $_[4];
    my $use_sim_tree      = $_[5];
    my $use_probabilistic = $_[6];
    my $mloc_Opts         = $_[7];
    my $dpDir             = $_[8];
    my $resDir            = $_[9];

    my $new_only = 1;    ## compute only new aligs

    mkdir($tgtDir);

    open( TREE, "$tree_file" );
    my $tree_str = <TREE>;
    close(TREE);
    chomp($tree_str);

    ## get sim-tree structure, $tree is hash-ref, keys are nodeIDs in postorder from 0
    ## nodes: postordered list of node-labels (leaf-names & '$$nodesym')

    my ( $nodes, $tree ) = newick_tree_to_postorder2($tree_str);
    my ( $childs2, $nodes2 ) = getNodeLeafs( $tree, 0 );
    my $subtrees_simtree_href = getSubtrees( $tree, $min_seqs, $max_seqs )
      ;    ## key == value for max subtrees

    my %imMap = ();
    my $imIdx = 1;
    foreach my $idx ( 1 .. @{$nodes2} ) {
        next if $tree->{ $nodes2->[ $idx - 1 ] }->{LEAF};
        $imMap{ $nodes2->[ $idx - 1 ] } = $imIdx;
        $imIdx++;

    }

    print "\nmap node-id (preorder) -> intermediate-id\n";
    foreach my $key ( sort keys %imMap ) {
        print "imMap $key -> " . $imMap{$key} . "\n";
    }
    mkdir("$tgtDir/log");
    mkdir($resDir);
    my @fa_tree = read_fasta_file($fa_tree_file);
    #
    foreach my $id ( sort { $a <=> $b } keys %{$subtrees_simtree_href} ) {

        ## continue only for maximal trees
        next if ( $subtrees_simtree_href->{$id} != $id );

        print "\ntree $id used " . $subtrees_simtree_href->{$id} . "\n";

        my $subtree_dir = "$tgtDir/$id";
        mkdir("$tgtDir/$id");

        my ( $childs, $nodes ) = getNodeLeafs( $tree, $id );

        my $subset_fa =
          writeSubsetFasta( \@fa_tree, $childs, "$tgtDir/$id/$id.fa", 1 );

        ####don't need this so far
        #   if ($use_sim_tree) {
        #     print "use sim-tree structure for subtree node $id\n";
        #     my $subtree_str = subtree2Newick( $tree, $id, 1 );
        #     print "sim-tree for subtree node $id: $subtree_str\n";
        #     open( TREE, ">$tgtDir/$id/$id.subtree" );
        #     print TREE $subtree_str . ";";
        #     close(TREE);
        #   }
        #
        ## reuse dot plots, does not work locarnaP currently!
        if ( !$use_probabilistic ) {
            my $loc_pp_dir = "$subtree_dir/input";
            mkdir($loc_pp_dir);
            foreach my $key ( @{$childs} ) {
                system("ln -f  $dpDir/$key $loc_pp_dir/$key")
                  if ( -e "$dpDir/$key" );
            }
        }

        my $call = "mlocarna --tgtdir $subtree_dir ";
        $call .= "--treefile $tgtDir/$id/$id.subtree " if ($use_sim_tree);
        $call .= "$mloc_Opts --verbose ";

        if ($use_probabilistic) {

            $call .= "--probabilistic " if ( $call !~ /--probabilistic/ );
            system( $call
                  . " $tgtDir/$id/$id.fa 1>$tgtDir/log/$id.locarna_call 2>$tgtDir/log/$id.locarna_call_err"
            ) if ( !-e "$tgtDir/$id/results/result.aln" || !$new_only );
        }
        else {

            system( $call . " $tgtDir/$id/$id.fa" )
              if ( !-e "$tgtDir/$id/results/result.aln" || !$new_only );

        }

        print "eval subtree $id from simtree\n";
        open( TREE, "$tgtDir/$id/results/result.tree" );
        my $eval_tree_str = <TREE>;
        chomp($eval_tree_str);
        my ( $nodeList, $evalTree ) = newick_tree_to_postorder2($eval_tree_str);
        my ( $leafs_subtree, $nodes_subtree ) = getNodeLeafs( $evalTree, 0 );

        print "here eval subtree $id from simtree nodes "
          . join( ":", @{$nodes} ) . "\n";

        my $iidx = 0;

        ## go through nodes of subtree in postorder
        foreach my $node ( @{$nodes_subtree} ) {

            ## ignore leafs
            next if ( $evalTree->{$node}->{LEAF} );
            $iidx++ if ( !$evalTree->{$node}->{LEAF} );

            ## ignore subtrees out of range, basically only too small ones
            next
              if ( $evalTree->{$node}->{NLEAFS} < $min_seqs
                || $evalTree->{$node}->{NLEAFS} > $max_seqs );

            ## analyse subtree
            #      my ( $leafsH, $nodesH ) = getNodeLeafs( $evalTree, $node );
            my ( $childs3, $nodes3 ) = getNodeLeafs( $tree, $node + $id );
            print "eval node $node ("
              . $id
              . ") intermediate $iidx orig-node:"
              . ( $node + $id )
              . " leafs "
              . join( ":", @{$childs3} ) . "\n";

            ## intermediate name: either original intermediate idx if same subtree structure is used
            ## or if new structure is used: "intermediate_rootIntermID_idx.aln"
            my $im_filename;
            if ($use_sim_tree) {
                $im_filename = "intermediate_" . $imMap{ $node + $id } . ".aln";
            }
            else {
                $im_filename = "intermediate_$id\_$iidx.aln";
            }
            ##print "node=$node id=$id iidx:$iidx  im name: $im_filename\n";
            system(
"cp $tgtDir/$id/intermediates/intermediate$iidx.aln $resDir/$im_filename"
            );

            ## write out parent leafs for matrix file
            open( PL, ">$resDir/$im_filename.subtree_ids" );
            print PL join( " ", sort @{$childs} ) . "\n";
            close(PL);

            if ( -e "$tgtDir/$id/results/result.matrix" ) {
                system(
"cp $tgtDir/$id/results/result.matrix $resDir/$im_filename.matrix"
                );
            }

            ##ali fold
            my $tmp_dir = "$tgtDir/alifold_$$";
            my $currDir = getcwd;
            mkdir($tmp_dir);
            chdir($tmp_dir);
            my @call_alifold = readpipe(
"RNAalifold -r --noLP --color --aln ../../$tgtDir/$id/intermediates/intermediate$iidx.aln 2>/dev/null"
            );
            my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
            chomp($aln_cons_struct);
            $aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
            $aln_cons_struct = $1;
            open( CONS, ">../../$resDir/$im_filename.alifold" );
            print CONS $call_alifold[0];
            print CONS "$aln_cons_struct\n";

            print CONS $2;
            close(CONS);

            system("mv alirna.ps ../../$resDir/$im_filename.alirna.ps");
            system("mv aln.ps ../../$resDir/$im_filename.ps");
            chdir($currDir);
            system("rm -R $tmp_dir");
        }

    }
}

########################################ä
## newick_tree_to_postorder($tree string)
##
## Translates a newick tree format string into a list of nodes in postorder
##
## $tree string in newick tree format
##
## returns list of nodes/leaves in postorder, use $node_sym for inner nodes
##
## dies if tree is not parsable
##
########################################
sub newick_tree_to_postorder2 {
    my ($tree) = @_;

    my @list;

    #$tree =~ s/:[\d.e-]+//g;

    $tree =~ s/;$//;    # allow that the tree string is terminated by ';'

    my $brcount = 0;

    my $preID = -1;
    my @last  = ();
    my %tree  = ();

    for ( my $i = 0 ; $i < length $tree ; $i++ ) {
        my $c = substr $tree, $i, 1;

        if ( $c eq "(" ) {

            $brcount++;
            $preID++;
            push( @last, $preID );

            # print "open preID $preID\n";

            $tree{$preID}->{LEAF}   = 0;
            $tree{$preID}->{NLEAFS} = 0;
            $tree{$preID}->{CHILDS} = [];
            $tree{$preID}->{NAMES}  = [];

        }
        elsif ( $c eq ")" ) {

            my $l = pop(@last);

            if ($l) {

                #print "close preID $preID last $l\n";
                push( @{ $tree{ $last[$#last] }->{CHILDS} }, $l );
                push(
                    @{ $tree{ $last[$#last] }->{NAMES} },
                    @{ $tree{$l}->{NAMES} }
                );
                $tree{ $last[$#last] }->{NLEAFS} += @{ $tree{$l}->{NAMES} };
            }

            ## parse branch length of inner node
            my $rest = substr $tree, $i + 1;

            my $len = 0;
            my $t   = "";
            if ( $rest =~ /^([^(),]+)/ ) {
                $t = $1;
                $t =~ /\:(.*)/;
                $len = $1;
            }
            $tree{$l}->{BRANCH} = $len;
            $i += ( length $t );

            $brcount--;

            if ( $brcount < 0 ) {
                die "Parse error in tree.";
            }

            push @list, $node_sym;
        }
        elsif ( $c eq "," ) {
            ## ignore, although we could do syntax checking
        }
        else {
            my $rest = substr $tree, $i;
            $rest =~ /^([^(),]+)/;
            my $t = $1;

            my @part  = split( ":", $t );
            my $token = $part[0];
            my $len   = $part[1];

            if ( !$token ) {
                $t =~ /\:(.*)/;
                $len = $1;
            }

            #print "t $t $len last $last[$#last]\n";

            if ($token) {

                $preID++;

                $tree{$preID}->{LEAF}   = 1;
                $tree{$preID}->{NAME}   = $token;
                $tree{$preID}->{BRANCH} = $len;

                #print "token $token preID $preID last " . $last[$#last] . "\n";
                push( @{ $tree{ $last[$#last] }->{CHILDS} }, $preID );
                push( @{ $tree{ $last[$#last] }->{NAMES} },  $token );
                $tree{ $last[$#last] }->{NLEAFS}++;
            }

            #print "part $t token<".$token."> len<$len>\n";

            $i += ( length $t ) - 1;
            push @list, $token if ($token);
        }
    }

    # printTree( \%tree );
    return ( \@list, \%tree );
}

sub getNodeLeafs {
    my $tree   = $_[0];
    my $nodeID = $_[1];

    my @leafs;    ## leafs in postorder
    my @nodes;    ## postorder of all internal/leaf nodes

    ## inner nodes
    if ( !$tree->{$nodeID}->{LEAF} ) {
        my @childs = @{ $tree->{$nodeID}->{CHILDS} };

        my ( $childLeafs, $childNodes ) = getNodeLeafs( $tree, $childs[0] );
        push( @leafs, @{$childLeafs} );
        push( @nodes, @{$childNodes} );

        ( $childLeafs, $childNodes ) = getNodeLeafs( $tree, $childs[1] );
        push( @leafs, @{$childLeafs} );
        push( @nodes, @{$childNodes} );
        push( @nodes, $nodeID );          # if (!$tree->{$child}->{LEAF});
    }

    ## leafs
    if ( $tree->{$nodeID}->{LEAF} ) {
        return ( [ $tree->{$nodeID}->{NAME} ], [$nodeID] );
    }

    return ( \@leafs, \@nodes );
}

sub getSubtrees {
    my $tree     = $_[0];
    my $min_seqs = $_[1];
    my $max_seqs = $_[2];

    my %used_IDs = ();

    foreach my $nodeID ( sort { $a <=> $b } keys %{$tree} ) {

        my ( $childs1, $nodes1 ) = getNodeLeafs( $tree, $nodeID );

#  print "check node $nodeID leaf:" . $tree->{$nodeID}->{LEAF} . " nleafs:" . $tree->{$nodeID}->{NLEAFS} . " childs:" . join( ":", @{$childs1} ) . "\n";
        next if ( exists $used_IDs{$nodeID} );

        next
          if ( $tree->{$nodeID}->{LEAF}
            || $tree->{$nodeID}->{NLEAFS} < $min_seqs
            || $tree->{$nodeID}->{NLEAFS} > $max_seqs );

        print "use Node $nodeID\n";
        my ( $childs, $nodes ) = getNodeLeafs( $tree, $nodeID );

        $used_IDs{$nodeID} = $nodeID;

        foreach my $node ( @{$nodes} ) {
            $used_IDs{$node} = $nodeID
              if ( !$tree->{$node}->{LEAF}
                && !exists $used_IDs{$node}
                && $tree->{$node}->{NLEAFS} >= $min_seqs
                && $tree->{$node}->{NLEAFS} <= $max_seqs );
        }

        print "node $nodeID leafs "
          . join( ":", @{$childs} )
          . " subnodes "
          . join( ":", @{$nodes} ) . "\n";

    }

    return \%used_IDs;
}

sub writeSubsetFasta {
    my $fa_aref  = $_[0];
    my $ids_aref = $_[1];
    my $out_file = $_[2];
    my $info     = $_[3];

    #  my $seqs_href   = $_[0];
    #  my $out_file    = $_[1];
    #  my $header_href = $_[2];

    open( OUT, ">$out_file" ) or die "Cannot open file $out_file! Exit...\n\n";

    ## order is important !!!
    foreach my $key ( sort { $a <=> $b } @{$ids_aref} ) {

        die "Error! Seq-ID $key does not esist in fasta! Exit...\n\n"
          if ( !exists $fa_aref->[0]->{$key} );

        print OUT ">$key " . $fa_aref->[2]->{$key} . "\n";
        print OUT $fa_aref->[0]->{$key} . "\n";

        next if ( !$info );
        ## write meta info: given structure
        print OUT $fa_aref->[3]->{$key}->{"#FS"} . " #FS\n"
          if ( exists $fa_aref->[3]->{$key}->{"#FS"} );
        print OUT $fa_aref->[3]->{$key}->{"#S"} . " #S\n"
          if ( exists $fa_aref->[3]->{$key}->{"#S"} );
    }

    ## order is important !!!
    #  foreach my $key ( sort { $a <=> $b } keys %{$seqs_href} ) {
    #    print OUT ">$key";
    #    print OUT " " . $header_href->{$key} if ($header_href);
    #    print OUT "\n";
    #
    #    print OUT $seqs_href->{$key} . "\n";
    #  }

    close(OUT);
    return $out_file;
}

sub subtree2Newick {
    my $tree   = $_[0];
    my $nodeID = $_[1];
    my $root   = $_[2];

    my $str = "";

    if ( !$tree->{$nodeID}->{LEAF} ) {
        my @childs = @{ $tree->{$nodeID}->{CHILDS} };

        if ( $tree->{ $childs[0] }->{LEAF} ) {
            $str .= $tree->{ $childs[0] }->{NAME};
        }
        else {
            $str .= "(" . subtree2Newick( $tree, $childs[0] ) . ")";
        }

        if ( $tree->{ $childs[1] }->{LEAF} ) {
            $str .= "," . $tree->{ $childs[1] }->{NAME};
        }
        else {
            $str .= ",(" . subtree2Newick( $tree, $childs[1] ) . ")";
        }
    }
    $str = "(" . $str . ")" if ($root);

    return $str;
}

sub rankSubtrees {
    my $subtreeDir = $_[0];
    my $treeDir    = $_[1];
    my $evaluate   = $_[2];

    my @subtrees = ();

    open( TREE, "$treeDir" );
    my $tree_str = <TREE>;
    close(TREE);
    chomp($tree_str);

    ## get sim-tree structure, $tree is hash-ref, keys are nodeIDs in postorder from 0
    my ( $nodes, $tree ) = newick_tree_to_postorder2($tree_str);
    my ( $childs2, $nodes2 ) = getNodeLeafs( $tree, 0 );
    foreach my $file ( glob("$subtreeDir/*.aln") ) {

        my %node = ();
        my $aln  = read_clustalw_alnloh($file);
        $node{FILE}      = $file;
        $node{BRANCHSUM} = 0;
        $node{BRANCH}    = 0;
        $node{DENSITY}   = 0;
        $node{NAMES}     = [ map { $_->{name} } @{$aln} ];
        $node{NODEID}    = 0;
        $node{NLEAFS}    = @{$aln};


        ## fill node with scores
        NodeScore( \%node, $aln, $childs2, $tree_matrix );

        if ($evaluate) {
            ( $node{CLASS}, $node{QUAL}, $node{QUAL_ABS} ) = evalNode( \%node );
        }
        else {
            $node{CLASS}    = $node{NAMES};
            $node{QUAL}     = "0.0";
            $node{QUAL_ABS} = "0.0";
        }

        push( @subtrees, \%node );

    }
    return \@subtrees;
}

sub read_clustalw_alnloh {
    my ($aln_filename) = @_;

    my $fh;

    open( $fh, "$aln_filename" )
      || die "Can not read alignment file $aln_filename\n";
    my ( $aln, $names, $header ) = read_clustalw_alignment($fh);
    close $fh;

    if ( !$header ) { die "Missing CLUSTAL header in $aln_filename."; }

    ## generate list of hash alignment
    my @alnloh;

    for my $name (@$names) {
        my $seq = $aln->{$name};
        my $entry = { name => $name, seq => $seq };
        push @alnloh, $entry;
    }

    return \@alnloh;
}

sub makeCol {
    my $col_file = $_[0];

    #  my $tmp_dir = "_col_tmp_" . $$;
    #  mkdir($tmp_dir);
    #  system("column -t $col_file > $tmp_dir/t");
    #  system("mv $tmp_dir/t $col_file");
    #  system("rm -R -f $tmp_dir");

    ## read content and spit into columns
    open( my $in_file, $col_file );
    my @cont    = ();
    my $max_col = 0;
    while ( my $line = <$in_file> ) {
        chomp($line);
        my @t = split( " ", $line );
        $max_col = scalar(@t) if ( $max_col < scalar(@t) );
        push( @cont, \@t );
    }
    close($in_file);

    ## get max length of all columns
    my @len = ();
    map { push( @len, 0 ) } 1 .. $max_col;

    foreach my $ci ( 0 .. $max_col - 1 ) {
        foreach my $l (@cont) {
            if ( @{$l} >= $ci + 1 ) {
                my $clen = length( $l->[$ci] );
                $len[$ci] = $clen if ( $len[$ci] < $clen );
            }
        }
    }

    ## print out cols with fixed width
    my $col_spc = 2;

    open( my $out, ">$col_file.col_$$" );

    foreach my $l (@cont) {
        map { printf $out ( "%-*s", $len[$_] + $col_spc, $l->[$_] ) }
          0 .. @{$l} - 1;
        print $out "\n";
    }
    close($out);

    system("mv $col_file.col_$$ $col_file");
}

sub read_clustalw_alignment {
    my ($fh) = @_;

    my %aln;

    my @names = ();    ## keep a list of names for order

    my $line;

    my $clustal_header = 0;
    if ( ( $line = <$fh> ) =~ /^CLUSTAL/ ) {
        $clustal_header = 1;
        $line           = <$fh>;
    }

    do {
        if ( $line =~ /^([^\s]+)\s+(.+)/ ) {
            my $name = $1;
            my $seq  = $2;

            if ( !exists $aln{$name} ) {
                push @names, $name;
            }

            $aln{$name} .= $seq;

        }
    } while ( $line = <$fh> );

    return ( \%aln, \@names, $clustal_header );
}

sub NodeScore {
    my $node             = $_[0];
    my $node_aln         = $_[1];
    my $tree_ids         = $_[2];
    my $tree_matrix_file = $_[3];

    $node->{SCI}    = 0;
    $node->{COV_EI} = 1.0;

    ## alfold consensus MFE
    $node->{MFE} = 0;
    if ( -e $node->{FILE} . ".alifold" ) {
        my @alifold = readpipe( "cat " . $node->{FILE} . ".alifold" );
        $alifold[2] =~ /\(\s*(\S+)\s*=/;
        $node->{MFE} = $1;
    }

    ## MPI
    my $mpi = 1;

    my @getmpi =
      readpipe("perl $path/scoreAln.pl -i $node->{FILE} -f CLUSTALW -s mpi");
    if (@getmpi) {
        $mpi = $getmpi[0];
        chomp($mpi) if ($mpi);
    }
    $node->{MPI} = $mpi;

    ## REL
    $node->{REL} = [ 1, 1, 1, 1 ];
    if ( -e $node->{FILE} . ".rel_signal" ) {

        open( IN, $node->{FILE} . ".rel_signal" );
        my @rel_dat = <IN>;
        close(IN);

        ## line $rel_dat[3] is: 'REL_SCORES 48.92:38.57:53.42:42.12' -> split 2nd filed into array-ref
        $node->{REL} = [ split( ":", ( split( " ", $rel_dat[3] ) )[1] ) ];
    }

    ## RNAz
    $node->{RNAZ} = [ 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1 ];
    my $use_rnaz = 1;
    if ($use_rnaz) {
        $node->{RNAZ} = getRNAzScores( $node->{FILE}, "RNAz", "-n" );

        ## overwrite alifold mfe with RNAz consensus mfe if we have strange values
        if ( $node->{MFE} >= 0 && $node->{RNAZ}->[4] < 0 ) {
            $node->{MFE} = $node->{RNAZ}->[4];
        }

        if (   $node->{RNAZ}->[3] < 0
            && $node->{MFE} < 0
            && $node->{RNAZ}->[6] < 0
            && $node->{RNAZ}->[5] < 0 )
        {
            $node->{SCI} = abs( $node->{MFE} / $node->{RNAZ}->[3] );
        }
        else {
            $node->{SCI} = 0.01;
            $node->{MFE} = 0.0;
        }

        if ( $node->{RNAZ}->[6] < 0 ) {
            ## new 0.4
            $node->{COV_EI} = ( abs( $node->{RNAZ}->[6] ) )**0.5;
            ## new 0.7
            $node->{COV_EI} = 1 if ( $node->{COV_EI} < 1 );
        }

    }

    ## score
    $node->{SCORE}     = 0;
    $node->{SCORE_OLD} = 0;
    $node->{CENTERC}   = 0;

    $node->{LEN_MEAN} = 0;
    $node->{LEN_SD}   = 0;
    $node->{LEN_ALN}  = length( $node_aln->[0]->{seq} );
    $node->{SVM}      = 0;

    ## OVERLAP
    $node->{OVERLAP} = check_node_overlap( $node->{NAMES}, "$data_map", 0.1 );

    ## Greylist
    $node->{GL_HIT} = 0;

    ## score simple is just average of score submatrix for current subtree
    my %leaf2idx = ();    ## stores score-matrix row for each leaf-id

    @{$tree_ids} = sort { $a <=> $b } @{$tree_ids};

    map { $leaf2idx{ $tree_ids->[$_] } = $_ } 0 .. $#{$tree_ids};
    my @leaf_idx =
      map { $leaf2idx{$_} if ( exists $leaf2idx{$_} ); } @{ $node->{NAMES} };

    my $sc = getMatrixSum( $tree_matrix_file, \@leaf_idx );
    $node->{SCORE} = $sc;
    print "\nmatrix avg = $sc " . ( $sc / ( $node->{NLEAFS} - 1 ) ) . "\n";
    print $node->{NAMES} . " " . $node->{FILE} . "\n";

    ## get locarna subtree score
    $node->{SCORE_LOC} = 1;
    if ( -e $node->{FILE} . ".matrix" ) {
        %leaf2idx = ();
        my $sub_ids = readSubset( $node->{FILE} . ".subtree_ids", 1 );
        @{$sub_ids} = sort { $a <=> $b } @{$sub_ids};
        map { $leaf2idx{ $sub_ids->[$_] } = $_ } 0 .. $#{$sub_ids};
        @leaf_idx = ();
        @leaf_idx = map { $leaf2idx{$_} if ( exists $leaf2idx{$_} ); }
          @{ $node->{NAMES} };
        my $sc_loc = getMatrixSum( $node->{FILE} . ".matrix", \@leaf_idx );

     #$node->{SCORE_LOC} =  ($sc_loc)/( $node->{NLEAFS}-1);
     #$node->{SCORE_LOC} =  $sc_loc/($center_subtree_max-@{ $node->{NAMES} }+1);

        if ( -e $node->{FILE} . ".rel_signal" ) {

            #$sc_loc = $sc_loc / ( $node->{NLEAFS} - 1 );
            $sc_loc = $sc_loc;
        }
        $node->{SCORE_LOC} = $sc_loc;
        print "score locarna matrix = $sc_loc \n";
    }

    #print "idx:".join(":",@leaf_idx)."\n";
    #print "matrix_score:".$sc."\n";

    ## FINAL SCORE
    if ( $center_tree_type == 3 ) {
        ## tree based on kernel sim matrix
        $sc = $sc / ( $node->{NLEAFS} - 1 );
    }
    else {
        $sc = $sc / ( $node->{NLEAFS} - 1 );
    }
    $node->{SIZE_BAL} = 1 / ( ( $node->{MPI} * 100 )**( 1 / $node->{NLEAFS} ) );
    ## full_05_new_04
    $node->{SCORE_SIMPLE} =
      $node->{SCORE_LOC} *
      $node->{SIZE_BAL} *
      $node->{SCI}**0.7 *
      $node->{COV_EI};
    ##new
#$node->{SCORE_SIMPLE} =  $node->{SCORE_LOC} * $node->{SIZE_BAL} * $node->{SCI} * $node->{SIZE_BAL} * $node->{COV_EI};

}

sub getRNAzScores {
    my $aln_file  = $_[0];
    my $rnaz_path = $_[1];
    my $rnaz_opts = $_[2];

    my @rnaz_scores = ( 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 );

    $rnaz_opts = "" if ( !$rnaz_opts );

    #die "Cannot find RNAz $rnaz_path! Exit...\n" if ( !-e $rnaz_path );
    my $call = readpipe("$rnaz_path $rnaz_opts $aln_file");

    if ( $call =~ /Mean pairwise identity:\s+(\S+)/ ) {
        $rnaz_scores[0] = $1;
        $call =~ /Shannon entropy:\s+(\S+)/;
        $rnaz_scores[1] = $1;
        $call =~ /G+C content:\s+(\S+)/;
        $rnaz_scores[2] = $1;
        $call =~ /Mean single sequence MFE:\s+(\S+)/;
        $rnaz_scores[3] = $1;
        $call =~ /Consensus MFE:\s+(\S+)/;
        $rnaz_scores[4] = $1;
        $call =~ /Energy contribution:\s+(\S+)/;
        $rnaz_scores[5] = $1;
        $call =~ /Covariance contribution:\s+(\S+)/;
        $rnaz_scores[6] = $1;
        $call =~ /Combinations\/Pair:\s+(\S+)/;
        $rnaz_scores[7] = $1;
        $call =~ /Mean z-score:\s+(\S+)/;
        $rnaz_scores[8] = $1;
        $call =~ /Structure conservation index:\s+(\S+)/;
        $rnaz_scores[9] = $1;
        $call =~ /SVM RNA-class probability:\s+(\S+)/;
        $rnaz_scores[10] = $1;
    }
    else {
        die "Error! Unexpected RNAz output:\n$call\n";
    }

    return \@rnaz_scores;
}

sub getMatrixSum {
    my $mat_file = $_[0];
    my $idx      = $_[1];
    my $currDir  = getcwd;

    open( IN, "$mat_file" );
    my @rows = <IN>;
    @rows =
      map { my @a = split( " ", $rows[$_] ); $rows[$_] = \@a } 0 .. $#rows;
    close(IN);

    @{$idx} = sort { $a <=> $b } @{$idx};

    my @row_sel = ();
    my $mode    = 1;
    if ( $mode == 1 ) {
        ## extract submatrix for all idx
        @row_sel =
          map { my $a = $rows[$_]; my @b = @{$a}[ @{$idx} ]; \@b } @{$idx};
    }
    else {
        ## extract full lines for all idx
        @row_sel = map { my $a = $rows[$_]; $a } @{$idx};
    }

    #print "sel idx:" . join( ":", @{$idx} ), "\n";
    my $sum_all = 0;
    my $num     = @{$idx};
    foreach my $id ( 0 .. $#{$idx} ) {

        print "row idx" . $id . " val:" . join( ":", @{ $row_sel[$id] } ) . " ";

        my $sum = 0;
        map { $sum += $_ } @{ $row_sel[$id] };

        print " sum $sum\n";
        $sum_all += $sum;
    }
    my $avg = $sum_all / ($num);

    #print "matrix avg $avg\n";
    return $avg;
}

sub readSubset {
    my $file     = $_[0];
    my $sub_idx  = $_[1];
    my $max_rank = $_[2];

    open( IN, "$file" ) or die "Cannot open file $file! Exit...\n\n";

    my @set = ();

    my $idx = 0;
    while ( my $line = <IN> ) {
        $idx++;
        next if ( $idx != $sub_idx );
        chomp $line;
        @set = split( " ", $line );
        last;
    }
    close(IN);

    print
"(readSubset) File $file contains less than $sub_idx lines! \n"
      if ( !@set );

    if (@set) {
        $max_rank = -1 if ( !$max_rank );
        $max_rank = @set if ( $max_rank == -1 || $max_rank > @set );

        @set = @set[ 0 .. ( $max_rank - 1 ) ];
    }
    return ( \@set );
}

sub check_node_overlap {
    my $node_ids    = $_[0]; ## leaf ids in one node
    my $id_map_file = $_[1]; ## usually data.map file line: "122 SEQ23#20#100#+"
    my $black_overlap =
      $_[2];    ## minimal required overlap of two frags to become blacklisted

    ## map id to frag
    my $frags   = read_fragments($id_map_file);
    my %id2frag = ();
    map { $id2frag{ $_->{VALUE} } = $_ } @{$frags};

    ## matrix frags
    my $node_frags = [];
    map { push( @{$node_frags}, $id2frag{$_} ) } @{$node_ids};

    ## matrix frag ols
    @{$node_frags} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$node_frags};
    my $matrix_ols =
      fragment_overlap( $node_frags, $node_frags, $black_overlap, 1 );

    my $node_overlap = 0;
    $node_overlap = 1 if ( @{$matrix_ols} );

    return $node_overlap;
}

sub read_fragments {
    my $frag_file = $_[0];

    my @frags = ();
    open( IN, $frag_file )
      or die "Cannot find fragment file $frag_file. Exit...\n\n";
    while ( my $line = <IN> ) {
        chomp $line;

        my @ent = split( " ", $line );
        die "Wrong fragment hash format! Exit...\n\n" if ( @ent < 2 );

        my $frag = str2frag( $ent[1] );
        $frag->{VALUE} = $ent[0];

        push( @frags, $frag );
    }
    close(IN);

    @frags = sort { $a->{SEQID} cmp $b->{SEQID} } @frags;

    return \@frags;
}

sub str2frag {
    my $str = $_[0];

    ## str: frag-key (eg. "SEQ1#2#13#+")
    $str =~ s/\s+//g;
    my @key = split( "#", $str );
    die
"Wrong fragment key format! see: $str - Expected 'SEQ1#5#15#+' Exit...\n\n"
      if ( @key != 4 );

    my $frag = {};
    $frag->{SEQID}  = $key[0];
    $frag->{START}  = $key[1];
    $frag->{STOP}   = $key[2];
    $frag->{STRAND} = $key[3];
    $frag->{KEY}    = $str;

    return $frag;
}

sub fragment_overlap {
    my $fh1           = $_[0];
    my $fh2           = $_[1];
    my $ol_cutoff     = $_[2];
    my $ignore_strand = $_[3];

    my @overlaps = ();

    @{$fh1} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$fh1};
    @{$fh2} = sort { $a->{SEQID} cmp $b->{SEQID} } @{$fh2};

    my $f2_curr_start_idx = 0;
    my $f2_curr_end_idx   = 0;
    my @f2_keys           = ();

 #    map { print "f1 order: $_ " . $fh1->[$_]->{SEQID} . "\n"; } 0 .. $#{$fh1};
 #    map { print "f2 order: $_ " . $fh2->[$_]->{SEQID} . "\n"; } 0 .. $#{$fh2};

    ## check overlaps in more efficient way than all-vs-all
    ## achive this by special sorting of frags
    foreach my $f1 ( 0 .. $#{$fh1} ) {

#     print "$f1 $f2_curr_start_idx f1:" . $fh1->[$f1]->{KEY} ." ". $fh2->[$f2_curr_start_idx]->{KEY}." $f2_curr_start_idx ".@{$fh2}."\n";

        ## check if current (implies all) fragments in fh2 are "greater" that current f1
        next
          if ( $f2_curr_start_idx >= @{$fh2}
            || $fh1->[$f1]->{SEQID} lt $fh2->[$f2_curr_start_idx]->{SEQID} );

        ## get all fragments on the same seq in set $fh2
        if ( $fh1->[$f1]->{SEQID} ge $fh2->[$f2_curr_start_idx]->{SEQID} ) {

            # print "range: $f2_curr_start_idx .. $f2_curr_end_idx\n";
            # my $oldend = $f2_curr_end_idx;

            ## set new start-index in $fh2 for correct SEQID
            ## stop at index with equal SEQID or "lower" SEQID (than $f1 SEQID is not in $fh2 )
            while ($f2_curr_start_idx < @{$fh2}
                && $fh1->[$f1]->{SEQID} ne $fh2->[$f2_curr_start_idx]->{SEQID}
                && $fh1->[$f1]->{SEQID} gt $fh2->[$f2_curr_start_idx]->{SEQID} )
            {
                $f2_curr_start_idx++;
            }

            $f2_curr_end_idx = $f2_curr_start_idx;
            if (   $f2_curr_start_idx >= @{$fh2}
                || $fh1->[$f1]->{SEQID} lt $fh2->[$f2_curr_start_idx]->{SEQID} )
            {
                next;
            }

            ## set new end-index in $fh2 for current SEQID
            while ($f2_curr_end_idx < @{$fh2} - 1
                && $fh1->[$f1]->{SEQID} eq
                $fh2->[ $f2_curr_end_idx + 1 ]->{SEQID} )
            {
                $f2_curr_end_idx++;
            }

            ## sort frags according to startpos
            @f2_keys = sort { $fh2->[$a]->{START} <=> $fh2->[$b]->{START} }
              $f2_curr_start_idx .. $f2_curr_end_idx;

#map {print "skip $_ ".$fh2->[$_]->{SEQID}."\n";} $oldend+1 .. $f2_curr_start_idx-1;
# print "range: $f2_curr_start_idx .. $f2_curr_end_idx\n";
        }

        ## check all pairs, sorting allows to stop if f2 frag starts after frag f1
        foreach my $f2 (@f2_keys) {

#     print "f1:" . $fh1->[$f1]->{SEQID} . " f2:" . $fh2->[$f2]->{SEQID} . "\n";

            next
              if ( $fh1->[$f1] == $fh2->[$f2] )
              ;    ## check if we compare the same object (pointer)

            ## done: if add_reverse=true then cmsearch_scan_reverse should be false, otherwise double hits
            ## wrong because we scon only on original input fasta,
            ## not on additional reverse seqs which were added as fragments if wanted
            ## check for same strand, be careful with cmsearch hits on reverse strand
            next
              if (!$ignore_strand
                && $fh1->[$f1]->{STRAND} ne $fh2->[$f2]->{STRAND} );

            my $start1 = $fh1->[$f1]->{START};
            my $start2 = $fh2->[$f2]->{START};
            my $end1   = $fh1->[$f1]->{STOP};
            my $end2   = $fh2->[$f2]->{STOP};

            next if ( $end2 < $start1 );    ## f2 frag ist before f1
            last
              if ( $end1 < $start2 )
              ; ## f2 is after f1, due to sorting no overlap is possible anymore

            my $len1 = $end1 - $start1 + 1;
            my $len2 = $end2 - $start2 + 1;

            # 1-2-2-1 ; 1-2-1-2; 2-1-2-1; 2-1-1-2
            my $left  = $start2 - $start1;
            my $right = $end1 - $end2;
            my $overlap_len =
              ( abs( $len1 + $len2 ) - abs($left) - abs($right) ) / 2;

            ## use shortest frag as reference length
            my $ol_ref_len = min( $len1, $len2 );
            my $overlap_ratio = $overlap_len / $ol_ref_len;
            $overlap_ratio = sprintf( "%.2f", $overlap_ratio );

            next if ( $overlap_ratio < $ol_cutoff );

       #      print "$f1 $f2  $fh1->[$f1]->{KEY} $fh2->[$f2]->{KEY} overlap!\n";
            push( @overlaps, [ $f1, $f2, $overlap_ratio ] );
        }
    }
    return \@overlaps;
}

sub getNodeInfo {
    my $node = $_[0];
    my $nidx = $_[1];
    my $rank = $_[2];

    my $info_str = "";

    $info_str .= "cent $nidx rank $rank";
    $info_str .= " score " . sprintf( "%.3f", $node->{SCORE} );
    $info_str .= " cent " . $node->{CENTERC};

    $info_str .= " SVM " . $node->{SVM};

    $info_str .= " class " . join( ":", @{ $node->{CLASS} } );
    $info_str .= " names " . join( ":", @{ $node->{NAMES} } );

    $info_str .= " sci " . sprintf( "%.2f", $node->{SCI} );
    $info_str .= " MPI " . $node->{MPI};
    $info_str .= " REL " . join( ":", @{ $node->{REL} } );
    $info_str .= " mfe " . $node->{MFE};
    $info_str .= " dens " . $node->{DENSITY};
    $info_str .= " BSum " . $node->{BRANCHSUM};
    $info_str .= " nodeID " . $node->{NODEID};
    $info_str .= " nleafs " . $node->{NLEAFS};

    my $tmp = $1
      if ( $node->{FILE} =~ /intermediate_(\d+\.aln)/
        || $node->{FILE} =~ /intermediate_(\d+_\d+\.aln)/ );
    $info_str .= " file " . $1;
    $info_str .= " RNAZ " . join( ":", @{ $node->{RNAZ} } );
    $info_str .= " LEN_MEAN " . $node->{LEN_MEAN};
    $info_str .= " LEN_SD " . $node->{LEN_SD};
    $info_str .= " LEN_ALN " . $node->{LEN_ALN};
    $info_str .= " QUAL " . $node->{QUAL};
    $info_str .= " QUAL_ABS " . $node->{QUAL_ABS};
    $info_str .= " SC_SIMPLE " . sprintf( "%.2f", $node->{SCORE_SIMPLE} );
    $info_str .= " OL " . $node->{OVERLAP};
    $info_str .= " GL " . $node->{GL_HIT};
    $info_str .= " LOC " . sprintf( "%.2f", $node->{SCORE_LOC} );
    $info_str .= " COV_EI " . sprintf( "%.2f", $node->{COV_EI} );
    $info_str .= " SIZE_BAL " . sprintf( "%.2f", $node->{SIZE_BAL} );

    return $info_str;

}

sub aln2alifold {
    my $aln_file  = $_[0];
    my $tmp_path  = $_[1];
    my $vrna_path = $_[2];

    ## alifold result to get consensus structure string for infernal and some nice pictures
    my $tmp_dir = "$tmp_path/alifold_$$";
    my $currDir = getcwd;

    system("mkdir $tmp_dir");
    system("cd $tmp_dir");

    #mkdir($tmp_dir);
    chdir($tmp_dir);
    system("RNAalifold -r --noLP --color --aln $aln_file 2>/dev/null");

#my @call_alifold = readpipe( "RNAalifold -r --noLP --color --aln $aln_file 2>/dev/null" );
# my $aln_cons_struct = $call_alifold[1];    ## store cons-struct
# chomp($aln_cons_struct);
# $aln_cons_struct =~ /([\(\).]*)(.*)/;      ## extract cons-struct
# $aln_cons_struct = $1;
# open( CONS, ">$aln_file.alifold" );
# print CONS $call_alifold[0];
# print CONS "$aln_cons_struct\n";
# print CONS $2;
# close(CONS);
# system("mv alirna.ps $aln_file.alirna.ps");
# system("mv aln.ps $aln_file.ps");
    chdir($currDir);

    #system("rm -R $tmp_dir");
}
