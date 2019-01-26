# AnnotateRnaz.py ---
#
# Filename: AnnotateRnaz.py
# Description:
# Author: Joerg Fallmann
# Maintainer:
# Created: Sat Jan 26 12:45:25 2019 (+0100)
# Version:
# Package-Requires: ()
# Last-Updated: Sat Jan 26 16:21:54 2019 (+0100)
#           By: Joerg Fallmann
#     Update #: 185
# URL:
# Doc URL:
# Keywords:
# Compatibility:
#
#

# Commentary:
# This script is a replacement for rnazAnnotate.pl
# rnazAnnotate can not handle the output from version 2 adequatly
# This script uses the bedtools API to fast intersect an annotation Bed
# with output from RNAz

# Change Log:
#
#
#
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with GNU Emacs.  If not, see <http://www.gnu.org/licenses/>.
#
#

# Code:

#!/usr/bin/env python3

import sys
import glob
import argparse
from io import StringIO
import gzip
import traceback as tb
from pybedtools import *
import pybedtools
import re
import tempfile

def parseargs():
    parser = argparse.ArgumentParser(description='Intersect RNAz output with Annotation from BED')
    parser.add_argument("-b", "--bed", type=str, help='Annotation BED file')
    parser.add_argument("-i", "--input", type=str, help='RNAz output')
    parser.add_argument("-o", "--bedout", type=str, help='Annotation BED output')
    parser.add_argument("-r", "--rnazout", type=str, help='Annotation rnaz output')
    return parser.parse_args()

def annotate(bed, input, bedout, rnazout):
    try:

        pybedtools.set_tempdir('.')  # Make sure we do not write somewhere we are not supposed to
        anno = BedTool(bed)
        rnaz=readrnaz(input)
        tmpbed = pybedtools.BedTool(rnaztobed(rnaz), from_string=True)

        intersection = tmpbed.intersect(anno,wa=True,wb=True,s=True)  # intersect strand specific, keep all info on a and b files

        bedtornaz(intersection, rnaz, bedout, rnazout)

        return 1

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        print(''.join(tbe.format()),file=sys.stderr)

def readin(file):
    try:
        if '.gz' in file:
            f = gzip.open(file,'rt')
        else:
            f = open(file,'rt')
        return f

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        print(''.join(tbe.format()),file=sys.stderr)

def readrnaz(rnaz):
    try:
        toparse = readin(rnaz)
        tointersect = {}
        header = []
        for line in toparse:
            if '#' in line[0]:
                tointersect['header']=line.strip()
                line = re.sub('^#','',line)
                cont = line.strip().split('\t')
                foi = cont.index('seqID') # need to find which column contains seqID
                sf = cont.index('start') # need to find which column contains start
                ef = cont.index('end') # need to find which column contains end
                if 'strand' in cont:# need to find which column contains strand
                    df = cont.index('strand')
                else:
                    df = None
            else:
                content = line.strip().split('\t')
                newid=re.split('\.|\,|\s|\\|\/|\_', content[foi])[1]  # I can only hope that we have species.chromosome.whatever as annotation in aln or maf, if not this is hardly parseable
                if df:
                    longid = '_'.join([newid, content[sf], content[ef], 'RNAzresult', '0', content[df]])
                    tointersect[longid] = content
                else:
                    longid = '_'.join([newid, content[sf], content[ef], 'RNAzresult', '0', '+'])
                    tointersect[longid] = content
                    longid = '_'.join([newid, content[sf], content[ef], 'RNAzresult', '0', '-'])
                    tointersect[longid] = content

        return tointersect

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        print(''.join(tbe.format()),file=sys.stderr)


def rnaztobed(rnaz):
    try:
        tmpbed = []
        for key in rnaz:
            if key != 'header':
                tmpbed.append('\t'.join(key.split('_')))

        return '\n'.join(tmpbed)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        print(''.join(tbe.format()),file=sys.stderr)

def bedtornaz(bed, rnaz, bedout, rnazout):
    try:
        b = open(bedout,'w')
        r = open(rnazout,'w')

        annotatedbed=[]
        annotatedrnaz=[]
        annotatedrnaz.append(str.join('\t',[rnaz['header'],'Annotation']))
        for line in open(bed.fn):
            out = line.strip().split("\t")
            annotatedbed.append(str.join('\t',out[0:3]+out[9:10]+out[4:6]))
            key = str.join('_',out[0:6])
            annotatedrnaz.append(str.join('\t',rnaz[key]+out[9:10]))

        print(str.join('\n', annotatedbed),file=b)
        print(str.join('\n', annotatedrnaz),file=r)

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        print(''.join(tbe.format()),file=sys.stderr)


def closefile(file):
    try:
        file.close()

    except Exception as err:
        exc_type, exc_value, exc_tb = sys.exc_info()
        tbe = tb.TracebackException(
            exc_type, exc_value, exc_tb,
        )
        print(''.join(tbe.format()),file=sys.stderr)




####################
####    MAIN    ####
####################
if __name__ == '__main__':
    args=parseargs()
    annotate(args.bed, args.input, args.bedout, args.rnazout)
######################################################################
# AnnotateRnaz.py ends here
