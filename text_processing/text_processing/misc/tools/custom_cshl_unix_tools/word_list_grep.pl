#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Std;

sub parse_command_line();
sub load_word_list();
sub compile_regex(@);
sub usage();

my $word_list_file;
my $input_file ;
my $output_file;
my $find_complete_words ;
my $find_inverse; 
my $find_in_specific_column ;
my $find_case_insensitive ;
my $skip_first_line ;


##
## Program Start
##
usage() if @ARGV==0;
parse_command_line();

my @words = load_word_list();

my $regex = compile_regex(@words);

# Allow first line to pass without filtering?
if ( $skip_first_line ) {
	my $line = <$input_file>;
	print $output_file $line ;
}


##
## Main loop
##
while ( <$input_file> ) {
	my $target = $_;


	# If searching in a specific column (and not in the entire line)
	# extract the content of that one column
	if ( $find_in_specific_column ) {
		my @columns = split ;

		#not enough columns in this line - skip it
		next if ( @columns < $find_in_specific_column ) ;

		$target = $columns [ $find_in_specific_column - 1 ] ;
	}

	# Match ?
	if ( ($target =~ $regex) ^ ($find_inverse) ) {
		print $output_file $_ ;
	}
}

close $input_file;
close $output_file;

##
## Program end
##


sub parse_command_line()
{
	my %opts ;
	getopts('siwvc:o:', \%opts) or die "$0: Invalid option specified\n";

	die "$0: missing word-list file name\n" if (@ARGV==0); 

	$word_list_file = $ARGV[0];
	die "$0: Word-list file '$word_list_file' not found\n" unless -e $word_list_file ;

	$find_complete_words = ( exists $opts{w} ) ;
	$find_inverse = ( exists $opts{v} ) ;
	$find_case_insensitive = ( exists $opts{i} ) ;
	$skip_first_line = ( exists $opts{s} ) ;


	# Search in specific column ?
	if ( defined $opts{c} ) {
		$find_in_specific_column = $opts{c};

		die "$0: invalid column number ($find_in_specific_column).\n"
			unless $find_in_specific_column =~ /^\d+$/ ;
			
		die "$0: invalid column number ($find_in_specific_column).\n"
			if $find_in_specific_column <= 0; 
	}
	else {
		$find_in_specific_column = 0 ;
	}


	# Output File specified (instead of STDOUT) ?
	if ( defined $opts{o} ) {
		my $filename = $opts{o};
		open $output_file, ">$filename" or die "$0: Failed to create output file '$filename': $!\n" ;
	} else {
		$output_file = *STDOUT ;
	}



	# Input file Specified (instead of STDIN) ?
	if ( @ARGV>1 ) {
		my $filename = $ARGV[1];
		open $input_file, "<$filename" or die "$0: Failed to open input file '$filename': $!\n" ;
	} else {
		$input_file = *STDIN;
	}
}

sub load_word_list()
{
	open WORDLIST, "<$word_list_file" or die "$0: Failed to open word-list file '$word_list_file'\n" ;
	my @words ;
	while ( <WORDLIST> ) {
		chomp ;
		s/^\s+//;
		s/\s+$//;
		next if length==0;
		push @words,quotemeta $_;
	}
	close WORDLIST;

	die "$0: Error: word-list file '$word_list_file' is empty!\n" 
       		unless @words;

	return @words;	
}

sub compile_regex(@)
{
	my @words = @_;

	my $regex_string = join ( '|', @words ) ;
	if ( $find_complete_words ) {
		$regex_string = "\\b($regex_string)\\b"; 
	}
	my $regex;

	if ( $find_case_insensitive ) {
		$regex = qr/$regex_string/i ;
	} else {
		$regex = qr/$regex_string/;
	}

	return $regex;
}

sub usage()
{
print <<EOF;

Word-List Grep
Copyright (C) 2009 - by A. Gordon ( gordon at cshl dot edu )

Usage: $0 [-o OUTPUT] [-s] [-w] [-i] [-c N] [-v] WORD-LIST-FILE [INPUT-FILE]

   -s   - do not filter first line - always output the first line from the input file.
   -w   - search for complete words (not partial sub-strings).
   -i   - case insensitive search.
   -v   - inverse - output lines NOT matching the word list.
   -c N - check only column N, instead of entire line (line split by whitespace).
   -o OUT - specify output file (default = STDOUT).
   WORD-LIST-FILE - file containing one word per line. These will be used
          for the search. 
   INPUT-FILE - (optional) read from file (default = from STDIN).



EOF

	exit;
}
