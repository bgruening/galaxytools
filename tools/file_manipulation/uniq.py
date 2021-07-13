import subprocess
import sys

"""
    We only need that file because galaxy do not understand the -t $'\t' term.
    Otherwise that would be the right XML-only solution:
    sort -u
        $ignore_case
        $is_numeric
        -t $'\t'
        #if $adv_opts.adv_opts_selector=="advanced":
            -k$adv_opts.column_start,$adv_opts.column_end
        #end if
        -o $outfile
        $input
"""

if sys.argv[1].strip() != "false":
    ignore_case = sys.argv[1]
else:
    ignore_case = ""

if sys.argv[2].strip() != "false":
    is_numeric = sys.argv[2]
else:
    is_numeric = ""

try:
    col_start = sys.argv[3]
    col_end = sys.argv[4]
    com = "sort -u %s %s -t '	' -k%s,%s -o %s %s" % (
        is_numeric,
        ignore_case,
        col_start,
        col_end,
        sys.argv[5],
        sys.argv[6],
    )
except Exception:
    # no advanced options selected
    com = "sort -u %s %s -t '	' -o %s %s" % (
        is_numeric,
        ignore_case,
        sys.argv[3],
        sys.argv[4],
    )

subprocess.call(com, shell=True)
