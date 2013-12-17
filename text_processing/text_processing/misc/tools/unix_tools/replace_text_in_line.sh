#!/bin/sh

##
## Galaxy wrapper for simple sed find&replace command
##

FIND_PATTERN="$1"
REPLACE_PATTERN="$2"
INPUT="$3"
OUTPUT="$4"

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
echo "Replacing: $FIND_PATTERN"
echo "With: $REPLACE_PATTERN"

sed -r --sandbox "s/$FIND_PATTERN/$REPLACE_PATTERN/g" "$INPUT" > "$OUTPUT"
