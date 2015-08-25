#!/bin/sh
#
# IMAGE cheat run script
#
#
#
source /etc/bash.bashrc

# software path
# this is the path where the IMAGE path is
# Please change it accordingly 
VELPATH=~/Desktop/IMAGE_version2/
SSAHADIR=~/Desktop/IMAGE_version2/
WALKPATH=~/Desktop/IMAGE_version2/

PATH=$WALKPATH:$SSAHADIR:$VELPATH:$PATH
export PATH


#options to add
SOLEXANAME=76bp
DIRNAME=iteration_
NUMITERATIONS=4

# the aligner used
ALIGNER=bwa # or ssaha


# unzip the fastq files
gunzip $SOLEXANAME\_1.fastq.gz
gunzip $SOLEXANAME\_2.fastq.gz


# make some format 
# do not modify anything here
01_prepare_new_read_placed.pl read.placed.original
02_prepare_new_contigs_fa.pl read.placed.new contigs.fa.original
ln -s read.placed.new read.placed
ln -s contigs.fa.new contigs.fa


# IMAGE script
# you can insert new parameters here
000_walk_all_raw.pl -prefix $SOLEXANAME -dir_prefix $DIRNAME -iteration 1 -all_iteration $NUMITERATIONS -aligner_aligner $ALIGNER



# commands available to add
#
#compulsory parameters:-prefix 2540_5   (Solexa prefix lane)
# -iteration 1     (current iteration)
# -all_iteration 5 (number of iterations to be run after 1 iteration)
# -dir_prefix ite  (prefix for the iterations)
#
#optional parameters: 
# -kmer 31 (kmer option in velvet)
# -toignore 10     (exclude extension of less than n bases)
# -overhang 50  (to not extend if overhang has > 50bp) 
# -vel_ins_len 500 (insert length specified in velvet)
# -aligner_aligner bwa  (default aligner)
# -bwa_index is  (used for different parameters of bwa)



