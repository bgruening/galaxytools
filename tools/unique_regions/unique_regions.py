#!/usr/bin/env python

import os
import argparse

def remove_first_line(s):
    return s[s.find('\n')+1:s.rfind('\n')]

def build_cont(prev1, prev2, prev3):
    return "%s\t%i\t%i\n" % (prev1, prev2, prev3)

def overlap(infile):
    overlapping = ""
    unique = ""
    prev_chromo = ""
    prev_start = 0
    prev_end = 0

    #iterating the bed file, line by line
    for line in open(infile, 'r'):
        split_line = line.split("\t")
        chromo = str(split_line[0])
        start = int(split_line[1])
        end = int(split_line[2])

        #verifying whether the two successive regions are on the same chromosome,
        #if not, the comparison will be skipped
        if chromo != prev_chromo:
            prev_chromo = chromo
            unique += build_cont(prev_chromo, prev_start, prev_end)
        else:
            #seperating the unique regions
            if start >= prev_end:
                unique += build_cont(prev_chromo, prev_start, prev_end)
            else:
                #seperating the overlapping regions
                overlapping += build_cont(prev_chromo, prev_start, prev_end)
        prev_end = end
        prev_start = start

    #writing the last region
    unique += build_cont(prev_chromo, prev_start, prev_end)
    unique = remove_first_line(unique)
    return unique, overlapping

def unique_regions(args):
    l = 0
    #retrieving the two output groups (unique and overlapping regions)
    unique, overlapping = overlap(args.input)
    unover_file_name = "unique_regions_"+str(l)+".bed"
    unover_file_name = os.path.join(args.output, unover_file_name)
    #writing the unique regions to a seperate file
    with open(unover_file_name, 'w') as f1:
        f1.write(unique)
    #writing the overlapping regions to a seperate files
    with open("temp_overlapping_regions.bed", 'w') as f2:
        f2.write(overlapping)

    #if there were any overlapping regions, the process will be repeated
    #starting from the overlapping regions file.
    #the process will be repeated until there is no more overlapping regions
    while overlapping != "":
        l += 1
        args.input = "temp_overlapping_regions.bed"
        unique, overlapping = overlap(args.input)

        unover_file_name = "unique_regions_"+str(l)+".bed"
        unover_file_name = os.path.join(args.output, unover_file_name)
        with open(unover_file_name, 'w') as f1:
            f1.write(unique)
        with open("temp_overlapping_regions.bed", 'w') as f2:
            f2.write(overlapping)
    os.remove("temp_overlapping_regions.bed")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output', help="Path to the output files", default="output")
    args = parser.parse_args()
    unique_regions(args)
