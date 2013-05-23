#!usr/bin/env python

import os, sys
import subprocess

"""
    OSRA_DATA_FILES is set during the toolshed Installation
    If it is not set, use the standard configuration of OSRA. 
    That means we need to delete argument 4-7.
    That script is a hack, because we do not know the content of OSRA_DATA_FILES at xml evaluation time.

    osra -f $oformat $infile 
        -l \$OSRA_DATA_FILES/spelling.txt -a \$OSRA_DATA_FILES/superatom.txt
        > $outfile
"""

if not os.path.exists(sys.argv[7]):
    # OSRA_DATA_FILES path is not set or the spelling file is not existent
    sys.argv.pop(7) # superatom.txt path
    sys.argv.pop(6) # -a
    sys.argv.pop(5) # speling.txt path
    sys.argv.pop(4) # -l

sys.argv[0] = 'osra'
subprocess.call(sys.argv, stdout=sys.stdout)


