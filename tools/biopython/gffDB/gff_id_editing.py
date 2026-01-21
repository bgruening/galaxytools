#!/usr/bin/env python
"""
Description:
    Adds ID tags to every gff entry.

Usage:
    gff_id_editing.py -o <input_gff_file> -n <output_gff_file>

"""
import os
import argparse
import re
import sys

def main(oldGff, newGff):
    counter = {}
    # open the gff file
    with open(oldGff,'r') as o_file:
        with open(newGff,'w') as n_file:
            lineNum = 0
            for line in o_file:
                lineNum += 1
                if not re.match(r"\#",line):
                    line = line.strip()
                    lineList = line.split('\t')
                    # check if the gff file is consistent
                    if len(lineList) < 9:
                        print("Error: GFF file line " + str(lineNum) + " has less than 9 entries")
                        sys.exit()
                    seqname, source, feature, start, end, score, strand, frame, attribute = lineList
                    # check if we have an ID in the entry
                    attList = attribute.split(';')
                    noID = True
                    parent = 'noParNoID'
                    for att in attList:
                        if re.search(r"=",att):
                            tag, value = att.split('=')
                            if 'ID' == tag:
                                # update ID counter list
                                IDvalue = counter.get(value, 0)
                                counter[value] = IDvalue + 1
                                noID = False
                            elif ('parent' == tag) or ('Parent' == tag):
                                parent = value
                                # use '-' to concatenate multiple parent features
                                parent = parent.replace(',', '-')
                    # add ID if we have none
                    if noID:
                        # create ID base
                        newID = parent + ':' + feature
                        # check if ID already in use and update
                        IDvalue = counter.get(newID, 0)
                        counter[newID] = IDvalue + 1
                        # add counter to ID
                        newID = 'ID=' + newID + '_' + str(IDvalue + 1)
                        # add new ID to attribute list
                        attribute = newID + ';' + attribute
                    # creat new entry
                    n_file.write(seqname + "\t" + source + "\t" + feature + "\t" + start + "\t" + end + "\t" + score + "\t" + strand + "\t" + frame + "\t" + attribute + "\n")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'gff_source_editing.py', description='Adds ID tags to every gff entry.', prefix_chars='-+', epilog="")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--old_gff', '-o', dest='oldGff', required=True, help='input GFF file')
    parser.add_argument('--new_gff', '-n', dest='newGff', required=True, help='output GFF file')

    options = parser.parse_args()
    main(options.oldGff, options.newGff)

