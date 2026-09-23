#!/usr/bin/env perl
## A wrapper script to call chewTree
use strict;
use warnings;
use Cwd;
use English;
use File::Copy;
use File::Basename;


# Parse arguments
my ($input,
    $input_type,
    $phantcec_dm,
    $phantcec_tree) = @ARGV;

# Run program
my $abs_path = Cwd::abs_path($PROGRAM_NAME);
my $scriptdir = dirname($abs_path);
runChewTree();
exit(0);

# Run runMentaLiST, only the tree part
sub runChewTree {
	if ($input_type ne "distance") {
      # calc distance matrix from database
      my $result = system("python $scriptdir/scripts/mlst-hash-distance.py -i $input -o $phantcec_dm");
    } else
	{
	  copy($input, $phantcec_dm) or die "Copy failed: $!";
	}
    # calc tree from distance matrix
    system("python $scriptdir/scripts/mentalist_tree $phantcec_dm > $phantcec_tree");
    return 0;
}
