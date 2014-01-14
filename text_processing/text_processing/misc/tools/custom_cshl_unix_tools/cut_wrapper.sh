#!/bin/sh

##
## Galaxy wrapper for cut command.
##

##
## command line arguments:
##   complement flag (might be empty string)
##   what to cut (fields or characters)
##   cut list (e.g. 1,2,3,4)
##   input_file
##   output_file

COMPLEMENT="$1"
CUTWHAT="$2"
CUTLIST="$3"
INPUT="$4"
OUTPUT="$5"

if [ -z "$OUTPUT" ]; then
	echo "This script should be run from inside galaxy!" >&2
	exit 1
fi

if [ ! -r "$INPUT" ]; then
	echo "error: input file ($INPUT) not found!" >&2
	exit 1
fi

# Messages printed to STDOUT will be displayed in the "INFO" field in the galaxy dataset.
# This way the user can tell what was the command
if [ -z "$COMPLEMENT" ]; then
	echo -n "Extracting "
else
	echo "Deleting "
fi

case $CUTWHAT in
	-f)	echo -n "field(s) "
		;;

	-c)	echo -n "character(s) "
		;;
esac

echo "$CUTLIST"


glxy-cut $COMPLEMENT $CUTWHAT $CUTLIST < $INPUT > $OUTPUT

exit
