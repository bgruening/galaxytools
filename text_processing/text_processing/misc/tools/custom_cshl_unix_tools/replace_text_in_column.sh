#!/bin/sh

##
## Galaxy wrapper for simple awk find&replace command
##

FIND_PATTERN="$1"
REPLACE_PATTERN="$2"
COLUMN="$3"
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
echo "Replacing: $FIND_PATTERN"
echo "With: $REPLACE_PATTERN"
echo "In column: $COLUMN"

#adapt to awk's quirks - to pass an acutal backslash - two backslashes are required (just like in a C string)
REPLACE_PATTERN=${REPLACE_PATTERN//\\/\\\\}
glxy-gawk -v OFS="\t" --re-interval --sandbox "{ \$$COLUMN = gensub( /$FIND_PATTERN/, \"$REPLACE_PATTERN\", \"g\", \$$COLUMN ) ; print \$0 ; }" "$INPUT" > "$OUTPUT"
