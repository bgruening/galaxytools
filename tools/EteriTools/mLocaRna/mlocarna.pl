use List::Util qw/ min max /;


my ($center_fa_file, $tree_file, $p, $max_diff_am, $tau, $max_diff, $path) = @ARGV;


$mloc_opts = "-p $p --max-diff-am $max_diff_am --tau $tau  --max-diff $max_diff --alifold-consensus-dp";

#$num = $center_fa_file;
#$num =~ s/\D+\z//;
#print $str;
#print "\n";

print "tree file = $tree_file\n";

my $cmd_stk = "perl $path/mloc2stockholm.pl -file $tree_file";
$cmd_stk .= " -split_input yes -con_struct $tree_file.alifold";

#system("$cmd_stk");# == 0 or die " mloc2stockholm command was unable to run to completion:\n\n";


system("mlocarna $mloc_opts --treefile $tree_file $center_fa_file --tgtdir mloc") == 0 or die " mloc command was unable to run to completion:\n\n";


#my $ids_ext = readSubset( "$CLUSTER_DIR/$clus_idx.cluster/center.ids.ext", 1 );
#my @fa_ext_all = read_fasta_file("$CLUSTER_DIR/$clus_idx.cluster/center.fa.ext");
#my $fa_ext_merged = mergeFrags( \@fa_ext_all );
#writeSubsetFasta( $fa_ext_merged, $fa_ext_merged->[1], "$model_dir/cmfinder.fa", 1 );
#system("cp $tree_dir/mloc/results/result.aln  $model_dir/model.tree.stk ");
