#!/usr/bin/env python

import sys, os
import argparse
import shlex
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input', help='Input file name')
parser.add_argument('-o1','--output1', help='tabular output file')
parser.add_argument('-o2', '--output2', help='images output tarball')
parser.add_argument('-p', '--partitionFunction', action='store_true', help='partition function')
parser.add_argument('-m', '--mea', type=float, help='mean ensemble accuracy')
parser.add_argument('-s', '--parameters', help='arguments')
args=parser.parse_args()

myinput = open(args.input)

specialParameters=""
if args.partitionFunction:
    specialParameters += ' --partfunc=1'
if args.mea:
    specialParameters += ' --MEA=%s' % args.mea

args.parameters = specialParameters + args.parameters

print args.parameters

p = subprocess.check_output(shlex.split('RNAfold '+args.parameters), stdin=myinput)

# process output into a tab-seperated file
# when no partition function is calculated only the first four columns are added
# 1st tab: comment from the fasta file
# 2nd tab: sequence
# 3rd tab: dot-bracket
# 4th tab: energy
#-----pf or mea?----------
# 5th tab: dot-bracket
# 6th tab: energy
# 7th tab: dot-bracket
# 8th tab: energy
#-----mea?----------------
# 9th tab: dot-bracket
#10th tab: energy
# 9/11th tab: last string
lines=p.split('\n')
o=lines[0].replace('\t',' ')
if args.mea:
    for x in range(1, len(lines)):
        if x % 7 == 6:
            o+='\t'+lines[x]
        if x % 7 == 5:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('{').rstrip('}')
            o+='\t'+st
        if x % 7 == 4:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('{').rstrip('}')
            o+='\t'+st
        if x % 7 == 3:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('[').rstrip(']')
            o+='\t'+st
        if x % 7 == 2:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('(').rstrip(')')
            o+='\t'+st
        if x % 7 == 1:
            o+='\t'+lines[x]
        if x % 7 == 0:
            o+='\n'+lines[x].replace('\t',' ')
elif args.partitionFunction:
    for x in range(1, len(lines)):
        if x % 6 == 5:
            o+='\t'+lines[x]
        if x % 6 == 4:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('{').rstrip('}')
            o+='\t'+st
        if x % 6 == 3:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('[').rstrip(']')
            o+='\t'+st
        if x % 6 == 2:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('(').rstrip(')')
            o+='\t'+st
        if x % 6 == 1:
            o+='\t'+lines[x]
        if x % 6 == 0:
            o+='\n'+lines[x].replace('\t',' ')
else:
    for x in range(1, len(lines)):
        if x % 3 == 2:
            st=lines[x].split(' ',1)[1]
            st=st.strip().lstrip('(').rstrip(')')
            o+='\t'+st
        if x % 3 == 1:
            o+='\t'+lines[x]
        if x % 3 == 0:
            o+='\n'+lines[x].replace('\t',' ')
out=open(args.output1,'w')
out.write(o)
out.close()

# postscript images are collected in a tar ball
p = subprocess.check_output("tar -c *.ps",shell=True)
out2 = open(args.output2,'w')
out2.write(p)
out2.close()





