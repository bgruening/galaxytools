#!/usr/bin/env python

"""
    Remaps a gff file back to genomic intervals using the glimmer3 prediction table.
    Author: Bjoern Gruening
"""

import sys



def main():

    glimmer_file = sys.argv[1]
    gff_file = sys.argv[2]
    new_gff_file = open(sys.argv[3], 'w')


    glimmer_dict = dict()
    for line in open(glimmer_file):
        line = line.strip()
        if not line or line.startswith('#') or line.startswith('>'):
            continue
        cols = line.split()
        glimmer_dict[cols[0]] = {'start': int(cols[1]), 'end': int(cols[2]), 'strand': int(cols[3])}


    for line in open(gff_file):
        if not line or line.startswith('#'):
            new_gff_file.write(line)
            continue
        cols = line.split('\t')
        if cols[0] in glimmer_dict:
            gff_entry = glimmer_dict[cols[0]]
            if gff_entry['strand'] < 0:
                #start column in gff, we are on the negative strand so inverse it
                new_start = str(gff_entry['start'] + int(cols[4]))
                new_end = str(gff_entry['start'] + int(cols[3]))
                cols[6] = '-'
                cols[3] = new_start
                cols[4] = new_end
            else:
                #start column in gff
                new_start = str(gff_entry['start'] + int(cols[3]))
                new_end = str(gff_entry['start'] + int(cols[4]))
                cols[6] = '+'
                cols[3] = new_start
                cols[4] = new_end
        new_gff_file.write('\t'.join(cols))

    new_gff_file.close()

if __name__ == '__main__':
    main()

