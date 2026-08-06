#!/usr/bin/env perl

use strict;
use Getopt::Long;
use Bio::SeqIO;
use Pod::Usage;

my ($fasta_file, $tab_file, $coverage_co, $length_co, $repeat_co, $out_filtered, $out_repeats, $out_norepeats,$coverage_length_co, $summary_out, $filtered_repeats, $help);

GetOptions(
	'c|coverage-cutoff=s'	       => \$coverage_co,
	'l|length-cutoff=s'	       => \$length_co,
        'e|coverage-length-cutoff=s'   => \$coverage_length_co,
	'r|repeat_cutoff=s'	       => \$repeat_co, 
	'i|input=s'		       => \$fasta_file,
	't|tab=s'		       => \$tab_file,
	'f|filtered-out=s'	       => \$out_filtered,
	'o|output-repeats=s'	       => \$out_repeats,
	'u|output-norepeats=s'	       => \$out_norepeats,
        'n|filtered-repeats=s'         => \$filtered_repeats,
        's|summary=s'                  => \$summary_out,
	'h|help'		       => \$help
);

pod2usage(-verbose => 2) if ($help);
print "A fasta file is required. Please enter a fasta file using the -i flag.\n" if (!$fasta_file);
print "A spades tabs file is required. Please enter a tabs file using the -t flag\n" if (!$tab_file);
pod2usage(1) unless $fasta_file && $tab_file;

if (!$coverage_co)
{
   $coverage_co = 0.33;
}
if (!$length_co)
{
   $length_co = 1000;
}
if (!$coverage_length_co)
{
   $coverage_length_co = 5000;
}
if (!$repeat_co)
{
   $repeat_co = 1.75;
}
if (!$out_filtered)
{
   $out_filtered = "Discarded_sequences.fasta";
   print "Discarded sequences will be printed out to $out_filtered\n";
}
if (!$out_repeats)
{
   $out_repeats = "Filtered_sequences_with_repeats.fasta";
   print "Filtered sequences with repeats will be printed out to $out_repeats\n";
}
if (!$out_norepeats)
{
   $out_norepeats = "Filtered_sequences_no_repeats.fasta";
   print "Filtered sequences without repeats will be printed out to $out_norepeats\n";
}
if (!$filtered_repeats)
{
   $filtered_repeats = "Repeat_sequences.fasta";
   print "Repeat sequences will be printed out to $filtered_repeats\n";
}

die ("No tab file specified") unless ($tab_file);
die ("No fasta file specified") unless ($fasta_file);

##Read tab file and discard rows with comments
open TAB, '<', $tab_file or die "Could not open tab file: $?";
open SEQIN, '<', $fasta_file or die "Could not open tab file: $?";
open SEQOUT_REP, '>', $out_repeats or die "Could not open file for writing: $?";
open SEQOUT_NOREP, '>', $out_norepeats or die "Could not open file for writing: $?";
open SEQOUT_FILT, '>', $out_filtered if ($out_filtered);
open SEQOUT_FILT_REP, '>', $filtered_repeats or die "Could not open file for writing: $?";
open SUMMARY, '>', $summary_out if ($summary_out);


my $avg_coverage = 0;
my $num_contigs = 0;
my $cutoff_coverage;
my $cutoff_repeats;
my @stats;


while (<TAB>)
{
	chomp;
	push @stats, $_ unless (/^#/);
}

#Calculate average coverage.
foreach my $stat(@stats)
{
	my ($length, $coverage);
	(undef,$length, $coverage) = split(/\t+/, $stat);
        die "length or coverage not defined at $stat\n" unless ($length && ($coverage ne '' && $coverage >= 0));
	if ($length >= $coverage_length_co)
	{
		$avg_coverage = $avg_coverage + $coverage;
		$num_contigs++;
	}
}

$avg_coverage = $avg_coverage / $num_contigs;
$cutoff_coverage = $avg_coverage * $coverage_co;
$cutoff_repeats = $avg_coverage * $repeat_co;

print SUMMARY "Filter SPAdes repeats Results Summary\n======================================\n\n" if ($summary_out);
print SUMMARY "Paramaters used:\nLength cutoff for calcularing average cutoff: $coverage_length_co\nCoverage cutoff ratio: $coverage_co\nRepeat cutoff ratio: $repeat_co\nLength cutoff: $length_co\n\n" if ($summary_out);

print SUMMARY "Calculations:\nAverage coverage: $avg_coverage\nCoverage cutoff: $cutoff_coverage\nRepeat cutoff: $cutoff_repeats\n\nFile headers:\n" if ($summary_out);

my ($header, $seq_id, $seq); 
my $repeated = 0;
my $valid = 0;

#Summary strings:
my $discarded = "";
my $repeats = "";
my $filtered_rep = "";
my $filtered_norep = "";

while (my $line = <SEQIN>)
{
	if ($line =~ />/)
	{
		chomp $line;
                #Get the sequence name to compare against tab file
		$header = $line;
		$seq_id = $line =~ /(\w+)_length/;
		$seq = "";
	
		my $stat = shift @stats;
		die "Less rows in tab than sequences in seq file" unless $stat;
		my($name, $length, $coverage) = split(/\t+/, $stat);
		die "name or length not defined at $stat\n" unless ($name && $length);
                die "coverage is not defined at $stat\n" unless ($coverage ne '' && $coverage >= 0);
		die "Unmatched names $header and $name\n" unless ($header =~ /$name/i);

               #Entry passes the length and coverage cutoffs?
		if ($length >= $length_co && $coverage >= $cutoff_coverage)
		{
			$valid = 1;
			#Repeats
			if ($coverage >= $cutoff_repeats)
			{
				my $num_repeats = int($coverage/$avg_coverage);
				$header = $header."(".$num_repeats." copies)";
				print SEQOUT_REP $header,"\n";
                                $filtered_rep = $filtered_rep.$header."\n";
                                print SEQOUT_FILT_REP $header, "\n";
                                $repeats = $repeats.$header."\n";
				$repeated = 1;
			}
			else
			{
				print SEQOUT_REP $header, "\n";
                                $filtered_rep = $filtered_rep.$header."\n";
				print SEQOUT_NOREP $header, "\n";
                                $filtered_norep = $filtered_norep.$header."\n";
				$repeated = 0;
			}
		}
		elsif ($out_filtered)
		{
			$valid = 0;
			print SEQOUT_FILT $header,"\n";
                        $discarded = $discarded.$header."\n";
		}
	}
	else
	{
		if ($valid)
		{
			print SEQOUT_REP $line;
			if (!$repeated)
			{
				print SEQOUT_NOREP $line;
			}
                        else
                        {
                              print SEQOUT_FILT_REP $line;
                        }
		}
		elsif ($out_filtered)
		{
			print SEQOUT_FILT $line;
		}
	}
	
}

close TAB;
close SEQIN;
close SEQOUT_REP;
close SEQOUT_NOREP;
close SEQOUT_FILT;
close SEQOUT_FILT_REP;


#Get summary info:
if ($summary_out)
{
   print SUMMARY "Filtered sequences (with repeats):\n$filtered_rep\n";
   print SUMMARY "Filtered sequences (no repeats):\n$filtered_norep\n";
   print SUMMARY "Repeat sequences:\n$repeats\n";
   if ($out_filtered)
   {
      print SUMMARY "Discarded sequences:\n$discarded\n"; 
   }

   close SUMMARY;
}

die "More rows in stats file than sequences in the fasta file\n" if (scalar(@stats) > 0);
exit 0;


__END__



=head1 NAME

	filter_spades_repeats.pl - Filters contigs or scaffolds based on contig length and detects contigs/scaffolds with very high coverage.



=head1 USAGE

	filter_spades_output.pl -i <contigs/scaffolds input>
                                -t <stats input>
                                -o <output fasta with repeats>
                                -u <output fasta without repeats>
                                
                                Optional:
                                -c <coverage cutoff ratio> (default 0.33) 
				-l <length cutoff> (default: 1000)
                                -e <length cutoff for average coverage calculation> (default: 5000)
				-r <repeat cutoff ratio> (default (1.75)
                                -n <filtered repeated sequences>
                                -f <discarded sequences>
                                -s <output summary file>

                                For more information:
                                -h


=head1 INPUT

=over 8

=item B<-i>B<--input>

Contigs/Scaffolds fasta file.

=item B<-t>B<--tab>

The tabular output file from SPAdes. This file should have the following format:

      #name length   coverage

      NODE_1   31438 24.5116
      
      NODE_2   31354 2316.96

      NODE_3   26948 82.3294

=item B<-o>B<--output-repeats>

Output fasta file including the contigs marked as repeated.

=item B<-u>B<--output-norepeats>

Output fasta file excluding the contigs marked as repeated.

=item B<-c>B<--coverage-cutoff>

Mininum coverage ratio. 

	coverage_theshold = average_coverage * minimum_coverage_ratio.

Any contigs/scaffolds with coverage below the coverage_theshold will be eliminated.

=item B<-l>B<--length-cutoff>

Mininum length. Contigs below this length will be eliminated.

=item B<-e>B<--coverage-length-cutoff>

Minimum length to use for average coverage calculations.

=item B<-r>B<--repeat-cutoff>

Minimum repeats ratio. 

	repeat_threshold = average_coverage * repeat_ratio. 

Any contigs with coverage below this threshold will be considered to be repeated


=item B<-f>B<--filtered-out>

If specified, filtered out sequences will be written to this file.

=item B<-s>B<--summary>

A summary of results

=back
=cut
