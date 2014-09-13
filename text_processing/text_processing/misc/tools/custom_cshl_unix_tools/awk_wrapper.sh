#!/bin/sh

##
## Galaxy wrapper for AWK command
##

##
## command line arguments:
##   input_file
##   output_file
##   awk-program
##   input-field-separator
##   output-field-separator

INPUT="$1"
OUTPUT="$2"
PROG="$3"

shift 3

if [ -z "$PROG" ]; then
	echo usage: $0 INPUTFILE OUTPUTFILE PROGRAM-FILE >&2
	exit 1
fi

if [ ! -r "$INPUT" ]; then
	echo "error: input file ($INPUT) not found!" >&2
	exit 1
fi

if [ ! -r "$PROG" ]; then
	echo "error: awk-program file ($PROG) not found!" >&2
	exit 1
fi

# Messages printed to STDOUT will be displayed in the "INFO" field in the galaxy dataset.
# This way the user can tell what was the command
cat "$PROG"

glxy-gawk --sandbox -v OFS="	" -v FS="	" --re-interval -f "$PROG" "$INPUT" > "$OUTPUT"
