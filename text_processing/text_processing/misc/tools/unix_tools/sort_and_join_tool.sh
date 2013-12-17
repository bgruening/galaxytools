#!/bin/sh
#
# Sort-and-Join wrapper
#
# This is a wrapper for Galaxy's sort-and-join (aka Easy-Join) tool.
#
# It can join two unsorted files, with all of GNU's join options, with header support
#
# All parameters must be supplied.
# the sort_and_join_tool.xml file takes care of that.


## For quicker operation, allow sort to use alot of memory
##
SORT_MEMORY_LIMIT="-S 2G"

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

##
## First, sort the two files, potentially with a header
##
TMP_SORT1=$(mktemp -t glx_sort_and_join_file_1.XXXXXXXXXXX) || exit
cat "$INPUT1" | glxy-sort-with-header "$HAVE_HEADER" -t "	" "-k${COLUMN1},${COLUMN1}b" ${SORT_MEMORY_LIMIT} > "$TMP_SORT1" || exit

TMP_SORT2=$(mktemp -t glx_sort_and_join_file_2.XXXXXXXXXXX) || exit
cat "$INPUT2" | glxy-sort-with-header "$HAVE_HEADER" -t "	" "-k${COLUMN2},${COLUMN2}b" ${SORT_MEMORY_LIMIT} > "$TMP_SORT2" || exit

##
## Then, join the sorted files
##
HEADER_ARG=""
if [ "x$HAVE_HEADER" == "x1" ]; then
	HEADER_ARG="--header"
fi

echo join -e "$EMPTY_STRING" --auto-format $IGNORE_CASE $JOINTYPE -1 "$COLUMN1" -2 "$COLUMN2"
glxy-join $HEADER_ARG \
	-t "	" \
	-e "$EMPTY_STRING" --auto-format \
	$JOINTYPE -1 "$COLUMN1" -2 "$COLUMN2" \
	"$TMP_SORT1" "$TMP_SORT2" > "$OUTPUT"

## Save this for later
JOIN_EXIT_CODE=$?

## Remote the two temporary files
rm -f "$TMP_SORT1" "$TMP_SORT2"

exit ${JOIN_EXIT_CODE}
