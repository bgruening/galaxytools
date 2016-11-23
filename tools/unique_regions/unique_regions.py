#!/usr/bin/env python

import os
import argparse

def remove_first_line(s):
    return s[s.find('\n')+1:s.rfind('\n')]

def overlap(infile):
    overlapping = ""
    unique = ""
    prev_chromo = ""
    prev_start = 0
    prev_end = 0

    #reading the input file
    with open(infile, 'r') as to_check:
        lines = to_check.readlines()

    #iterating the bed file, line by line
    for line in lines:
        split_line = line[:-1].split("\t")
        chromo = str(split_line[0])
        start = int(split_line[1])
        end = int(split_line[2])

        #verifying whether the two successive regions are on the same chromosome,
        #if not, the comparison will be skipped
        if chromo != prev_chromo:
            prev_chromo = chromo
            unique_regions = "%s\t%i\t%i\n" % (prev_chromo, prev_start, prev_end)
            unique += unique_regions
            prev_end = end
            prev_start = start
        else:
            #seperating the unique regions
            if start >= prev_end:
                unique_regions = "%s\t%i\t%i\n" % (prev_chromo, prev_start, prev_end)
                unique += unique_regions
                prev_end = end
                prev_start = start
            else:
                #seperating the overlapping regions
                overlapping_regions = "%s\t%i\t%i\n" % (prev_chromo, prev_start, prev_end)
                overlapping += overlapping_regions
                prev_end = end
                prev_start = start

    unique = remove_first_line(unique)
    unique_regions = "\n%s\t%i\t%i\n" % (chromo, start, end)
    unique += unique_regions
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
    print unique

    #if there were any overlapping regions, the process will be repeated
    #starting from the overlapping regions file.
    #the process will be repeated until there is no more overlapping regions
    while os.stat("temp_overlapping_regions.bed").st_size > 0:
        l += 1
        args.input = "temp_overlapping_regions.bed"
        unique, overlapping = overlap(args.input)

        unover_file_name = "unique_regions_"+str(l)+".bed"
        unover_file_name = os.path.join(args.output, unover_file_name)
        with open(unover_file_name, 'w') as f1:
            f1.write(unique)
        with open("temp_overlapping_regions.bed", 'w') as f2:
            f2.write(overlapping)
        print unique
    os.remove("temp_overlapping_regions.bed")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input')
    parser.add_argument('-o', '--output', help="Path to the output files", default="output")
    args = parser.parse_args()
    unique_regions(args)
