#!/bin/bash

##
## Galaxy wrapper for GREP command.
##

set -o pipefail

##
## command line arguments:
##   input_file
##   output_file
##   regex
##   COLOR or NOCOLOR
##   [other parameters passed on to grep]

INPUT="$1"
OUTPUT="$2"
REGEX="$3"
COLOR="$4"

shift 4

if [ -z "$COLOR" ]; then
	echo usage: $0 INPUTFILE OUTPUTFILE REGEX COLOR\|NOCOLOR [other grep patameters] >&2
	exit 1
fi

if [ ! -r "$INPUT" ]; then
	echo "error: input file ($INPUT) not found!" >&2
	exit 1
fi

# Messages printed to STDOUT will be displayed in the "INFO" field in the galaxy dataset.
# This way the user can tell what was the command
echo "grep" "$@" "$REGEX"

if [ "$COLOR" = "COLOR" ]; then
	#
	# What the heck is going on here???
	# 1. "GREP_COLORS" is an environment variable, telling GREP which ANSI colors to use.
	# 2. "--colors=always" tells grep to actually use colors (according to the GREP_COLORS variable)
	# 3. first sed command translates the ANSI color to a <FONT> tag with blue color (and a <B> tag, too)
	# 4. second sed command translates the no-color ANSI command to a </FONT> tag (and a </B> tag, too)
	# 5. last seds adds <html> tags

	GREP_COLORS="ms=31" grep --color=always -P "$@" -- "$REGEX" "$INPUT" | \
		grep -v "^\[36m\[K--\[m\[K$" | \
		sed -r 's/\[[0123456789;]+m\[K?/<font color="blue"><b>/g' | \
		sed -r 's/\[m\[K?/<\/b><\/font>/g' | \
		sed -r '1i<html><body><pre>' | sed -r '$a</pre></body></html>' > "$OUTPUT" || exit 1

elif [ "$COLOR" = "NOCOLOR" ]; then
	grep -P "$@" -- "$REGEX" "$INPUT" | grep -v "^--$" > "$OUTPUT" || exit 1
else
	echo Error: third parameter must be "COLOR" or "NOCOLOR" >&2
	exit 1
fi
