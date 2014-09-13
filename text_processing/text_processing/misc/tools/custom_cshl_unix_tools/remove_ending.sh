#!/bin/sh

# Version 0.1 ,  15aug08
# Written by Assaf Gordon (gordon@cshl.edu)
#

LINES="$1"
INFILE="$2"
OUTFILE="$3"

if [ "$LINES" == "" ]; then
	cat >&2 <<EOF 
Remove Ending Lines

Usage: $0 LINES [INFILE] [OUTFILE]

   LINES - number of lines to remove from the end of the file
   [INFILE] - input file (if not specified - defaults to STDIN)
   [OUTFILE]- output file (if not specified - defaults to STDOUT)

Input Example:

#Chr	Start	End
chr1	10	15
chr1	40	20
chr1	21	14
total   3 chromosomes

Removing 1 line (the last line) produces:

#Chr	Start	End
chr1	10	15
chr1	20	40
chr	14	21

Usage Example:
   
   \$ $0 1 < my_input_file.txt > my_output_file.txt

EOF
	
	exit 1
fi

#Validate line argument - remove non-digits characters
LINES=${LINES//[^[:digit:]]/}

#Make sure the line strings isn't empty
#(after the regex above, they will either contains digits or be empty)
if [ -z "$LINES" ]; then
	echo "Error: bad line value (must be numeric)" >&2
	exit 1
fi

# Use default (stdin/out) values if infile / outfile not specified
[ -z "$INFILE" ] && INFILE="/dev/stdin"
[ -z "$OUTFILE" ] && OUTFILE="/dev/stdout"

#Make sure the input file (if specified) exists.
if [ ! -r "$INFILE" ]; then
	echo "Error: input file ($INFILE) not found!" >&2
	exit 1
fi


# The "gunzip -f" trick allows
# piping a file (gzip or plain text, real file name or "/dev/stdin") to head
gunzip -f <"$INFILE" | head -n "-${LINES}" > "$OUTFILE"

