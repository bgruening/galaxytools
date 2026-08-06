#!/usr/bin/env python
"""
Description:
    Converts the source data of a gff file into a specified string.

Usage:
    gff_source_editing.py -o <input_gff_file> -n <output_gff_file> -s <string> 

"""
import os
import argparse
import re
import sys

def main(oldGff, newGff, sourceString, addString):
    # check source string
    if len( re.findall(r"\s", sourceString) ) > 0:
        print("Error: Source string contains whitespace characters or similar.")
        sys.exit()
    # open the gff file
    with open(oldGff,'r') as o_file:
        with open(newGff,'w') as n_file:
            lineNum = 1
            for line in o_file:
                if not re.match(r"\#",line):
                    lineList = line.split('\t')
                    # check if the gff file is consistent
                    if len(lineList) < 9:
                        print("Error: GFF file line " + str(lineNum) + " has less than 9 entries")
                        sys.exit()
                    seqname, source, feature, start, end, score, strand, frame, attribute = lineList
                    # create a new source or add to an existing source
                    if addString:
                        source += str(sourceString)
                    else:
                        source = str(sourceString)
                    n_file.write(seqname + "\t" + source + "\t" + feature + "\t" + start + "\t" + end + "\t" + score + "\t" + strand + "\t" + frame + "\t" + attribute)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'gff_source_editing.py', description='Converts the source data of a gff file into a specified string.', prefix_chars='-+', epilog="")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--old_gff', '-o', dest='oldGff', required=True, help='input GFF file')
    parser.add_argument('--new_gff', '-n', dest='newGff', required=True, help='output GFF file')
    parser.add_argument('--string', '-s', dest='sourceString', required=True, help='source string')
    parser.add_argument('--add', '-a', dest='addString', action='store_true', help='only add the string to the source')

    options = parser.parse_args()
    main(options.oldGff, options.newGff, options.sourceString, options.addString)

