#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Parse STAR Chimeric.out.junction and convert it to fusion junction file
that is acceptable for CIRCexplorer. This script was modified from
https://github.com/alexdobin/STAR/blob/master/extras/scripts/filterCirc.awk
"""

import sys
from collections import defaultdict


def main():
    if len(sys.argv) != 3:
        sys.exit('star_parse.py junc out')
    junc = defaultdict(int)
    with open(sys.argv[1], 'r') as junc_f:
        for line in junc_f:
            flag = int(line.split()[6])
            if flag < 0:
                continue
            chr1, site1, strand1, chr2, site2, strand2 = line.split()[:6]
            if chr1 != chr2 or strand1 != strand2:
                continue
            if strand1 == '+':
                start = int(site2)
                end = int(site1) - 1
            else:
                start = int(site1)
                end = int(site2) - 1
            if start > end:
                continue
            junc_id = '%s\t%d\t%d' % (chr1, start, end)
            junc[junc_id] += 1
    with open(sys.argv[2], 'w') as outf:
        for i, j in enumerate(junc):
            outf.write('%s\tFUSIONJUNC_%d/%d\t0\t+\n' % (j, i, junc[j]))


if __name__ == '__main__':
    main()
