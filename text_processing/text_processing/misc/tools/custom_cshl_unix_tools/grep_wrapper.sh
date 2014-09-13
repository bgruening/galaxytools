#!/bin/sh

##
## Galaxy wrapper for GREP command.
##

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

if [ "$COLOR" == "COLOR" ]; then
	#
	# What the heck is going on here???
	# 1. "GREP_COLORS" is an environment variable, telling GREP which ANSI colors to use.
	# 2. "--colors=always" tells grep to actually use colors (according to the GREP_COLORS variable)
	# 3. first sed command translates the ANSI color to a <FONT> tag with blue color (and a <B> tag, too)
	# 4. second sed command translates the no-color ANSI command to a </FONT> tag (and a </B> tag, too)
	# 5. htmlize_pre scripts takes a text input and wraps it in <HTML><BODY><PRE> tags, making it a fixed-font HTML file.

	GREP_COLORS="ms=31" grep --color=always -P "$@" -- "$REGEX" "$INPUT" | \
		grep -v "^\[36m\[K--\[m\[K$" | \
		sed -r 's/\[[0123456789;]+m\[K?/<font color="blue"><b>/g' | \
		sed -r 's/\[m\[K?/<\/b><\/font>/g' | \
		htmlize_pre.sh > "$OUTPUT"


	if (( $? ));  then exit; fi

elif [ "$COLOR" == "NOCOLOR" ]; then
	grep -P "$@" -- "$REGEX" "$INPUT" | grep -v "^--$" > "$OUTPUT"
	if (( $? ));  then exit; fi
else
	echo Error: third parameter must be "COLOR" or "NOCOLOR" >&2
	exit 1
fi

exit 0
