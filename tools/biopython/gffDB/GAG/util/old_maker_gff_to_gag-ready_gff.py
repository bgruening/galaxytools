#! /usr/bin/env python

# Intended to deal with one lonely file -- the .gff input for our Bdorsalis WGS submission
# This script replaces the ID field with the Name field and updates the Parent fields accordingly

# example input:
# scaffold00080	maker	gene	106151	109853	.	+	.	ID=1;Name=BDOR_007864
# scaffold00080	maker	mRNA	106151	109853	.	+	.	ID=2;Name=BDOR_007864-RA;Parent=1

# example output:
# scaffold00080	maker	gene	106151	109853	.	+	.	ID=BDOR_007864
# scaffold00080	maker	mRNA	106151	109853	.	+	.	ID=BDOR_007864-RA;Parent=BDOR_007864

import sys

def main():

    if len(sys.argv) != 2:
        sys.stderr.write("usage: old_maker_gff_to_gag-ready_gff.py <file.gff>\n")
        sys.exit()

    input_gff = sys.argv[1]
    id2name = {}

    with open(input_gff, 'r') as gff:
        for line in gff:
            if line.startswith('#'):
                continue
            splitline = line.split('\t')
            if len(splitline) != 9:
                sys.stderr.write("Oops, each line should have 9 fields. This line looks like: ")
                sys.stderr.write(line + "\nI'm giving up, sorry.\n")
                sys.exit()

            attributes = splitline[8]
            splitattr = attributes.split(';')
            parent = None
            id_ = splitattr[0].split('=')[1]
            name = splitattr[1].split('=')[1]
            if len(splitattr) > 2:
                parent = splitattr[2].split('=')[1].strip()
            id2name[id_] = name

            # Write updated line to stdout
            outline = ""
            for field in splitline[:8]:
                outline += field + '\t'
            newattributes = "ID=" + name
            if parent:
                newattributes += ";Parent=" + id2name[parent]
            outline += newattributes
            outline += '\n'
            sys.stdout.write(outline)




#################################################################################
if __name__ == '__main__':
    main()



