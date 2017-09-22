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

parser.add_argument("--raw_file", dest="raw_files", action='append',
                  help="The raw data files for the sequence.")

parser.add_argument("--processed_file", dest="processed_files", action='append',
                  help="The processed data files for the sequence.")

parser.add_argument("--title", dest="title",
                  help="The title for the sequence.")

parser.add_argument("--summary", dest="summary",
                  help="The summary for the sequence.")

parser.add_argument("--overall_design", dest="overall_design",
                  help="The overall design for the sequence.")

parser.add_argument("--contributor", dest="contributors", action='append',
                  help="The contributor for the sequence.")

parser.add_argument("--supplementary_file", dest="supplementary_file",
                  help="The supplementary file for the sequence.")

parser.add_argument("--sra_center_name_code", dest="sra_center_name_code",
                  help="The SRA Center Name Code for the sequence.")

parser.add_argument("--growth_protocol", dest="growth_protocol",
                  help="The growth protocol for the sequence.")

parser.add_argument("--treatment_protocol", dest="treatment_protocol",
                  help="The treatment protocol for the sequence.")

parser.add_argument("--extract_protocol", dest="extract_protocol",
                  help="The extract protocol for the sequence.")

parser.add_argument("--library_construction_protocol", dest="library_construction_protocol",
                  help="The library construction protocol for the sequence.")

parser.add_argument("--library_strategy", dest="library_strategy",
                  help="The library strategy for the sequence.")


parser.add_argument("--data_processing_step", dest="data_processing_steps", action='append',
                  help="A data processing step for the sequence.")

parser.add_argument("--genome_build", dest="genome_build",
                  help="The genome build for the sequence.")

parser.add_argument("--processed_data_files_format_content", dest="processed_data_files_format_content",
                  help="The processed data files format and content for the sequence.")


#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")
#parser.add_argument("--title", dest="title", default="",
#                  help="The title for the sequence..")

parser.add_argument("-o", "--output", dest="output_file", default="records.tsv",
                  help="The name of the output file.", metavar="file_name")

args = parser.parse_args()

with open(args.output_file, 'w') as tsvfile:
    writer = csv.writer(tsvfile, delimiter=',')
    writer.writerow(["SERIES"])
    writer.writerow(["title", args.title])
    writer.writerow(["summary", args.summary])
    writer.writerow(["overall design", args.overall_design])
    for contributor in args.contributors:
        writer.writerow(["contributor", contributor])
    writer.writerow(["supplementary file", args.supplementary_file])
    writer.writerow(["SRA_center_name_code", args.sra_center_name_code])
    writer.writerow([])
    writer.writerow(["SAMPLES"])
    writer.writerow(["Sample name", "title", "source name", "organism", "characteristics: cell type", "characteristics: passages", "characteristics: strain", "characteristics: ChIP antibody", "molecule", "description", "processed data file", "raw file", "raw file", "raw file", "raw file"])
    no_samples = 3
    for i in range(no_samples):
        writer.writerow([])
    writer.writerow([])
    writer.writerow(["PROTOCOLS"])
    writer.writerow(["growth protocol", args.growth_protocol])
    writer.writerow(["treatment protocol", args.treatment_protocol])
    writer.writerow(["extract protocol", args.extract_protocol])
    writer.writerow(["library construction protocol", args.library_construction_protocol])
    writer.writerow(["library strategy", args.library_strategy])
    writer.writerow([])
    writer.writerow(["DATA PROCESSING PIPELINE"])
    for data_processing_step in args.data_processing_steps:
        writer.writerow(["data processing step", data_processing_step])
    writer.writerow(["genome build", args.genome_build])
    writer.writerow(["processed data files format and content", args.processed_data_files_format_content])
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
tsvfile.close()