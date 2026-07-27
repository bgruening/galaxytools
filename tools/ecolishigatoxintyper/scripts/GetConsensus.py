#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import sys
import os
import subprocess
import shutil

def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input', dest='input', help='input alignment file')
    parser.add_argument('-o', '--output', dest='output', help='output consensus file')
    args = parser.parse_args()

    sequences = []
    varsequences = []
    # read input
    with open(args.input) as alignments:
        for alignment in alignments:
            if alignment[0] != ">":
                sequences.append(alignment.rstrip())
    numsequences = len(sequences)
    for j in range(0, numsequences + 1):
        varsequences.append("")
    lstnumvariants = []
    lstnumhyphens = []
    # loop over the columns
    for i in range(0, len(sequences[0])):
        variants = []
        numhyphens = 0
        # loop over the original rows to obtain variants
        for j in range(0, numsequences):
            if sequences[j][i:i+1] == "-":
                numhyphens = numhyphens + 1
            if sequences[j][i:i+1] not in variants and sequences[j][i:i+1] != "-":
                variants.append(sequences[j][i:i+1])
        lstnumvariants.append(len(variants))
        lstnumhyphens.append(numhyphens)
        # create varsequences with a template of the variants
        for j in range(0, numsequences):
            if lstnumhyphens[i] == 0:
                varsequences[j] = varsequences[j] + "-"
            elif lstnumvariants[i] < 2:
                varsequences[j] = varsequences[j] + "-"
            else:
                varsequences[j] = varsequences[j] + sequences[j][i:i+1]
        if lstnumvariants[i] == 1 and lstnumhyphens[i] > 0:
            varsequences[numsequences] = varsequences[numsequences] + variants[0]
        else:
            varsequences[numsequences] = varsequences[numsequences] + "-"
    # loop over the columns, apply single variant
    for i in range(0, len(sequences[0])):
        if lstnumvariants[i] == 1 and lstnumhyphens[i] > 0:
            # loop over all the rows to apply single variant to "-"
            for j in range(0, len(sequences)):
                if sequences[j][i:i+1] == "-":
                    lstsequence = list(sequences[j])
                    lstsequence[i] = varsequences[numsequences][i:i+1]
                    sequences[j] = ''.join(lstsequence)            
    # loop over the rows of the sequences, apply multiple variants
    for j in range(0, len(sequences)):
        # loop over the columns
        for i in range(0, len(sequences[0])):
            variants = []
            if sequences[j][i:i+1] == "-" and lstnumvariants[i] > 1:
                # loop over the rows of the varsequences
                for k in range(0, numsequences):
                    if varsequences[k][i:i+1] not in variants and varsequences[k][i:i+1] != "-":
                        variants.append(varsequences[k][i:i+1])
                        if len(variants) == 1:
                            lstsequence = list(sequences[j])
                            lstsequence[i] = variants[0]
                            sequences[j] = ''.join(lstsequence)
                        else:
                            lstsequence[i] = variants[len(variants) - 1]
                            sequences.append(''.join(lstsequence))
    # eliminate duplicate sequences
    sequences_unique = list(set(sequences))
    # write consensus sequences to output
    consensus = open(args.output, 'w')
    n = 1
    for sequence in sequences_unique:
        consensus.write(">consensus_" + str(n) + "\n")
        consensus.write(sequence + "\n")
        n = n + 1


if __name__ == "__main__":
    __main__()