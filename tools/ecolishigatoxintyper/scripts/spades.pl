#!/usr/bin/env perl
## A wrapper script to call spades.py and collect its output
use strict;
use warnings;
use File::Temp qw/ tempfile tempdir /;
use File::Copy;
use Getopt::Long;

# Parse arguments
my ($out_contigs_file,
    $out_contigs_stats,
    $out_scaffolds_file,
    $out_scaffolds_stats,
    $out_log_file,
    $new_name,
    @sysargs) = @ARGV;

## GetOptions not compatible with parsing the rest of the arguments in an array.
## Keeping the not-so-nice parse-in-one-go method, without named arguments.
# GetOptions(
#     'contigs-file=s'    => \$out_contigs_file,
#     'contigs-stats=s'   => \$out_contigs_stats,
#     'scaffolds-file=s'  => \$out_scaffolds_file,
#     'scaffolds-stats=s' => \$out_scaffolds_stats,
#     'out_log_file=s'    => \$out_log_file,
# );

# my @sysargs = @ARGV;

# Create temporary folder to store files, delete after use
#my $output_dir = tempdir( CLEANUP => 0 );
my $output_dir = 'output_dir';
# Link "dat" files as fastq, otherwise spades complains about file format

# Create log handle
open my $log, '>', $out_log_file or die "Cannot write to $out_log_file: $?\n";

# Run program
# To do: record time
runSpades(@sysargs);
collectOutput($new_name);
extractCoverageLength($out_contigs_file, $out_contigs_stats);
extractCoverageLength($out_scaffolds_file, $out_scaffolds_stats);
print $log "Done\n";
close $log;
exit 0;

# Run spades
sub runSpades {
    my $cmd = join(" ", @_) . " -o $output_dir";
    my $return_code = system($cmd);
    if ($return_code) {
	print $log "Failed with code $return_code\nCommand $cmd\nMessage: $?\n";
	die "Failed with code $return_code\nCommand $cmd\nMessage: $?\n";
    }
    return 0;
}

# Collect output
sub collectOutput{
    my ($new_name) = @_;
    
    # Collects output
    #if a new name is given for the contigs and scaffolds, change them before moving them
    if ( $new_name ne 'NODE') {
        renameContigs($new_name);
    }
    else {
        if ( -e "$output_dir/contigs.fasta") {
            move "$output_dir/contigs.fasta", $out_contigs_file;
        }
        if ( -e "$output_dir/scaffolds.fasta") {
            move "$output_dir/scaffolds.fasta", $out_scaffolds_file;
        }
    }

    open LOG, '<', "$output_dir/spades.log" 
	or die "Cannot open log file $output_dir/spades.log: $?";
    print $log $_ while (<LOG>);
    return 0;
}

#Change name in contig and scaffolds file
sub renameContigs{
    my ($name) = @_;

    open my $in, '<',"$output_dir/contigs.fasta" or die $!;
    open my $out,'>', $out_contigs_file;

    while ( my $line = <$in>) {
        #remove the NODE_ so we can rebuilt the display_id with our contig name with the contig number.
        #also move the remainder of the length
        if ( $line =~ />NODE_(\d+)_(.+)/) {
            $line = ">$name" . "_$1 $2\n";
        }
        print $out $line;
    }
    close $in;
    close $out;
    

    open $in, '<',"$output_dir/scaffolds.fasta" or die $!;
    open $out,'>', $out_scaffolds_file;

    while ( my $line = <$in>) {
        #remove the NODE_ so we can rebuilt the display_id with our contig name with the contig number.
        #also move the remainder of the length
        if ( $line =~ />NODE_(\d+)_(.+)/) {
            $line = ">$name" . "_$1 $2\n";
        }
        print $out $line;
    }
    close $in;
    close $out;

}


# Extract
sub extractCoverageLength{
    my ($in, $out) = @_;
    open FASTA, '<', $in or die $!;
    open TAB, '>', $out or die $!;
    print TAB "#name\tlength\tcoverage\n";
    while (<FASTA>){
	next unless /^>/;
	chomp;
	die "Not all elements found in $_\n" if (! m/^>(NODE|\S+)_(\d+)(?:_|\s)length_(\d+)_cov_(\d+\.*\d*)/);
	my ($name,$n, $l, $cov) = ($1,$2, $3, $4);
	print TAB "$name" . "_$n\t$l\t$cov\n";
    }
    close TAB;
}
