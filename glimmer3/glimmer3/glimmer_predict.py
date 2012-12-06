#!/usr/bin/env python
"""
Input: DNA Fasta File
Output: Tabular
Return Tabular File with predicted ORF's
Bjoern Gruening
"""
import sys, os
import tempfile
from random import Random
import string
import subprocess
import shutil

def __main__():
    
    genome_seq_file = sys.argv[1]
    outfile_path = sys.argv[2]
    outfile_ext_path = sys.argv[3]


    tag = ''.join(Random().sample(string.letters+string.digits, 12))
    tempdir = tempfile.gettempdir()

    #longorfs = os.path.join(tempdir, tag + ".longorf")
    trainingset = os.path.join(tempdir, tag + ".train")
    icm = os.path.join(tempdir, tag + ".icm")

    longorfs = tempfile.NamedTemporaryFile()
    trainingset = tempfile.NamedTemporaryFile()
    icm = tempfile.NamedTemporaryFile()


    #glimmeropts = "-o0 -g110 -t30 -l"
    glimmeropts = "-o%s -g%s -t%s" % (sys.argv[4], sys.argv[5], sys.argv[6])
    if sys.argv[7] == "true":
        glimmeropts += " -l"


    """
        1. Find long, non-overlapping orfs to use as a training set
    """
    subprocess.Popen(["tigr-glimmer", "long-orfs", "-n", "-t", "1.15",
        genome_seq_file, "-"], stdout = longorfs, 
        stderr = subprocess.PIPE).communicate()

    """
        2. Extract the training sequences from the genome file
    """
    subprocess.Popen(["tigr-glimmer", "extract", "-t",
        genome_seq_file, longorfs.name], stdout=trainingset, 
        stderr=subprocess.PIPE).communicate()

    """
        3. Build the icm from the training sequences
    """

    # the "-" parameter is used to redirect the output to stdout
    subprocess.Popen(["tigr-glimmer", "build-icm", "-r", "-"], 
        stdin=open(trainingset.name), stdout = icm, 
        stderr=subprocess.PIPE).communicate()

    """
        Run Glimmer3
    """
    b = subprocess.Popen(["tigr-glimmer", "glimmer3", glimmeropts, 
        genome_seq_file, icm.name, os.path.join(tempdir, tag)], 
        stdout = subprocess.PIPE, stderr=subprocess.PIPE).communicate()

    #shutil.copyfileobj
    shutil.copyfile( os.path.join(tempdir, tag + ".predict"), outfile_path )
    shutil.copyfile( os.path.join(tempdir, tag + ".detail"), outfile_ext_path )



if __name__ == "__main__" :
    __main__()
