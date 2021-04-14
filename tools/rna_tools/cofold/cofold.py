#!/usr/bin/env python

import sys, os
import argparse
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", help="Input file name")
parser.add_argument("-o1", "--output1", help="tabular output file")
parser.add_argument("-s", "--parameters", help="arguments")
args = parser.parse_args()

myinput = open(args.input)

parameters = args.parameters

# we assume that the param files are located next to the python dir
script_dir = os.path.dirname(os.path.realpath(__file__))
parameters = args.parameters.replace("-P ", "-P %s/" % script_dir)
p = subprocess.check_output(shlex.split("CoFold " + parameters), stdin=myinput)

lines = p.split("\n")
# FASTA header
o = lines[0].replace("\t", " ")

for x in range(1, len(lines)):
    if x % 3 == 2:
        [seq, st] = lines[x].split(" ", 1)
        st = st.strip().lstrip("(").rstrip(")")
        o += "\t" + seq + "\t" + st
    if x % 3 == 1:
        o += "\t" + lines[x]
    if x % 3 == 0:
        o += "\n" + lines[x].replace("\t", " ")
out = open(args.output1, "w")
out.write(o)
out.close()
