#!/bin/sh

#
# NOTE:
#  This is a wrapper for GNU's join under galaxy
#  not ment to be used from command line (if you're using the command line, simply run 'join' directly...)
#
# All parameters must be supplied.
# the join_tool.xml file takes care of that.

JOINTYPE="$1"
HAVE_HEADER="$2"
EMPTY_STRING="$3"
IGNORE_CASE="$4"

INPUT1="$5"
COLUMN1="$6"
INPUT2="$7"
COLUMN2="$8"
OUTPUT="$9"

if [ "$OUTPUT" == "" ]; then	
	echo "This script is part of galaxy. Don't run it manually.\n" >&2
	exit 1;
fi

echo join $OUTPUT_FORMAT -e "$EMPTY_STRING" $IGNORE_CASE $JOINTYPE -1 "$COLUMN1" -2 "$COLUMN2" 

HEADER_ARG=""
if [ "x$HAVE_HEADER" == "x1" ]; then
	HEADER_ARG="--header"
fi

glxy-join $HEADER_ARG \
	-t "	" \
	-e "$EMPTY_STRING" -o auto \
	$JOINTYPE -1 "$COLUMN1" -2 "$COLUMN2" \
	"$INPUT1" "$INPUT2" > "$OUTPUT"
