#!/usr/bin/env python
# -*- coding: UTF-8 -*-

__author__ = "Bjoern Gruening"
__version__ = "0.1"
__date__ = ""
__license__ = ''

"""
    These script converts a Blast XML Output to the NCBI Feature Table format,
    required for submission of Serquence data.
    http://www.ncbi.nlm.nih.gov/WebSub/index.cgi?tool=genbank
    TODO: include t-RNA results
"""
import os, sys
import argparse
import csv

parser = argparse.ArgumentParser()
parser.add_argument("-c", "--cutoff", dest="cutoff",
                  help="", metavar="CUTOFF", type=float)
parser.add_argument("-f", "--frequency", dest="fs",
                  help="", metavar="RATE", type=float)
parser.add_argument("-o", "--order", dest="order", default=5,
                  help="", metavar="ORDER", type=int)
parser.add_argument("-t", "--type", dest="filterType", default="LPF",
                  help=".", metavar="FILTER_TYPE")


args = parser.parse_args()

with open('records.tsv','w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter=',')
    writer.writerow(["SERIES"])
    writer.writerow(["title"])
    writer.writerow(["summary"])
    writer.writerow(["overall design"])
    writer.writerow(["contributor"])
    writer.writerow(["contributor"])
    writer.writerow(["supplementary file"])
    writer.writerow(["SRA_center_name_code"])
    writer.writerow([])
    writer.writerow(["SAMPLES"])
    writer.writerow(["Sample name", "title", "source name", "organism", "characteristics: cell type", "characteristics: passages", "characteristics: strain", "characteristics: ChIP antibody", "molecule", "description", "processed data file", "raw file", "raw file", "raw file", "raw file"])
    no_samples = 3
    for i in range(no_samples):
        writer.writerow([])
    writer.writerow([])
    writer.writerow(["PROTOCOLS"])
    writer.writerow(["growth protocol"])
    writer.writerow(["treatment protocol"])
    writer.writerow(["extract protocol"])
    writer.writerow(["library construction protocol"])
    writer.writerow(["library strategy"])
    writer.writerow([])
    writer.writerow(["DATA PROCESSING PIPELINE"])
    writer.writerow(["data processing step"])
    writer.writerow(["data processing step"])
    writer.writerow(["data processing step"])
    writer.writerow(["data processing step"])
    writer.writerow(["data processing step"])
    writer.writerow(["genome build"])
    writer.writerow(["processed data files format and content"])
    writer.writerow([])
    writer.writerow(["PROCESSED DATA FILES"])
    writer.writerow(["file name", "file type", "file checksum"])
    no_processed_data_files = 3
    for i in range(no_processed_data_files):
        writer.writerow([])
    writer.writerow(["RAW FILES"])
    writer.writerow(["file name", "file type", "file checksum", "instrument model", "read length", "single or paired-end"])
    no_raw_files = 8
    for i in range(no_raw_files):
        writer.writerow([])
    writer.writerow([])
    writer.writerow([])
    writer.writerow([])
    writer.writerow(["PAIRED-END EXPERIMENTS"])
    solid_experiments = False
    if (solid_experiments):
        writer.writerow(["file name 1", "file name 2", "file name 3", "file name 4", "average insert size", "standard deviation"])
    else:
        writer.writerow(["file name 1", "file name 2", "average insert size", "standard deviation"])
