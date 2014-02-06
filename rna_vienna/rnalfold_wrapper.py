# Filename: rnalfold_wrapper.py
# Author: Torsten Houwaart	
# Version: 29.01.2014
#
# This script accepts deals with the fact that the input for RNALfold comes from stdin
#
# -i		Input file

import sys
import os
import re
import string
import commands
from subprocess import Popen, PIPE
from tempfile import NamedTemporaryFile

# This function is exceedingly useful, perhaps package for reuse?
def getopts(argv):
    opts = {}
    while argv:
	if argv[0][0] == '-':
	    opts[argv[0]] = argv[1]
	    argv = argv[2:]
	else:
	    argv = argv[1:]
    return opts

def main():
    args = sys.argv[1:]

    try:
        opts = getopts(args)
    except IndexError:
        print "Usage:"
        print " -i sequence -o outputfile -L span -c converflag -g gquadflag -t temperature -d danglingmode -n lonepairflag -a noGUflag -b noClosingGUFlag"
        return 0
    
    sequence = opts.get("-i")
    if sequence == None:
        print "No sequence specified."
        return -1
    outputfile = opts.get("-o")
    if outputfile == None:
        print "No output file specified."
        return -2
    span = opts.get("-L")
    if span == None:
        print "No span size specified."
        return -3
    temperature = opts.get("-t")
    if temperature == None:
        print "No temperature specified."
        return -4
    noconvert = opts.get("-c")
    gquad = opts.get("-g")
    dangling = opts.get("-d")
    if dangling == None:
        print "No dangling parameter specified."
        return -5
  
    nlp = opts.get("-n")
    noGU = opts.get("-a")
    noClosingGU = opts.get("-b")

    flags = ""
    if noconvert == "true":
        flags = "--noconv"
    if gquad == "true":
        flags += " --gquad"
    if nlp == "true":
        flags += " --noLP"
    if noGU == "true":
        flags += " --noGU"
    if noClosingGU == "true":
        flags += " --noClosingGU"
       
    temperatures =" -T" + temperature
    danglings = " -d" + dangling
    spans = " -L" + span

    target = NamedTemporaryFile().name 
    open( target, 'w' ).write( sequence )

    #commandline = "cat %s > %s" % ( target, outputfile )
    commandline = "RNALfold %s %s %s %s < %s > %s" % ( temperatures, danglings, spans, flags, target, outputfile )
    
    # run command
    errorcode, stdout = commands.getstatusoutput(commandline)
    
    # return error code
    return errorcode

if __name__ == "__main__":
    main()
