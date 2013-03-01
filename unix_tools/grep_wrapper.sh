#!/bin/sh

##
## Galaxy wrapper for GREP command.
##

# The only truly portable way of getting a fully resolved path of a file.
SCRIPT=$(perl -MCwd -e '$f=Cwd::realpath($ARGV[0]); die unless -e $f; print $f' "$0") || exit 1
SCRIPT_DIR=$(dirname "$SCRIPT")
ANSI_SCRIPT="$SCRIPT_DIR/ansi2html.sh"
[ -e "$ANSI_SCRIPT" ] || { echo "Error: utility script ($ANSI_SCRIPT) not found." >&2 ; exit 1 ; }

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
	echo "usage: $0 INPUTFILE OUTPUTFILE REGEX COLOR|NOCOLOR [other grep patameters]" >&2
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
	GREP_COLOR='1;34' grep --color=always -P "$@" -- "$REGEX" "$INPUT" |
		sh "$ANSI_SCRIPT" > "$OUTPUT" || exit 1

elif [ "$COLOR" = "NOCOLOR" ]; then
	grep -P "$@" -- "$REGEX" "$INPUT" | grep -v "^--$" > "$OUTPUT" || exit 1
else
	echo Error: third parameter must be "COLOR" or "NOCOLOR" >&2
	exit 1
fi
