#!/bin/sh

##
## Galaxy wrapper for SED command
##

##
## command line arguments:
##   input_file
##   output_file
##   sed-program
##   [other parameters passed on to sed]

SILENT="$1"
INPUT="$2"
OUTPUT="$3"
PROG="$4"

if [ -z "$PROG" ]; then
	echo usage: $0 SILENT INPUTFILE OUTPUTFILE SED-PROGRAM-FILE >&2
	exit 1
fi

if [ ! -r "$INPUT" ]; then
	echo "error: input file ($INPUT) not found!" >&2
	exit 1
fi

if [ ! -r "$PROG" ]; then
	echo "error: sed-program file ($PROG) not found!" >&2
	exit 1
fi

SILENT_ARG=""
if [ "x$SILENT" = "x-n" ]; then
	SILENT_ARG="-n"
fi

# Messages printed to STDOUT will be displayed in the "INFO" field in the galaxy dataset.
# This way the user can tell what was the command
cat "$PROG"

sed -r --sandbox $SILENT_ARG -f "$PROG" "$INPUT" > "$OUTPUT"
