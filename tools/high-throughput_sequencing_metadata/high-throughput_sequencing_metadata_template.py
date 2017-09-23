#!/usr/bin/env python
# -*- coding: UTF-8 -*-

__author__ = "Hassan Gamaleldin"
__version__ = "0.1"
__date__ = ""
__license__ = ''

"""
    This tool produces a tsv file containing the High-throughput sequencing metadata template (version 2.1).
"""
import os, sys
import argparse
import csv
import hashlib

parser = argparse.ArgumentParser()


parser.add_argument("--processed_file", dest="processed_files", action='append',
                  help="The processed data files for the sequence.")

parser.add_argument("--processed_file_name", dest="processed_file_names", action='append',
                  help="The processed file name for the last processed file input.")

parser.add_argument("--raw_file", dest="raw_files", action='append',
                  help="The raw data files for the sequence.")

parser.add_argument("--raw_file_name", dest="raw_file_names", action='append',
                  help="The raw file name for the last raw file input.")

parser.add_argument("--raw_file_instrument_model", dest="raw_file_instrument_models", action='append',
                  help="The raw file instrument model for the last raw file input.")

parser.add_argument("--raw_file_read_length", dest="raw_file_read_lengths", action='append',
                  help="The raw file read length for the last raw file input.")

parser.add_argument("--raw_file_single_paired_choice", dest="raw_file_single_paired_choices", action='append',
                  help="The raw file single paired choice for the last raw file input.")


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


parser.add_argument("--sample_name", dest="samples_name", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_title", dest="samples_title", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_source_name", dest="samples_source_name", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_organism", dest="samples_organism", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_cell_type", dest="samples_cell_type", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_passages", dest="samples_passages", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_strain", dest="samples_strain", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_chip", dest="samples_chip", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_antibody", dest="samples_antibody", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_molecule", dest="samples_molecule", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_description", dest="samples_description", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_processed_input_files", dest="samples_processed_input_files", action='append',
                  help="Sample input item.")

parser.add_argument("--sample_raw_input_files", dest="samples_raw_input_files", action='append',
                  help="Sample input item.")



parser.add_argument("--paired_file_name1", dest="paired_file_names1", action='append',
                  help="First file name for paired-end experiments.")

parser.add_argument("--paired_file_name2", dest="paired_file_names2", action='append',
                  help="Second file name for paired-end experiments.")

parser.add_argument("--paired_file_name3", dest="paired_file_names3", action='append',
                  help="Thrid file name for paired-end experiments.")

parser.add_argument("--paired_file_name4", dest="paired_file_names4", action='append',
                  help="Fourth file name for paired-end experiments.")

parser.add_argument("--paired_avg", dest="paired_avgs", action='append',
                  help="Average insert size for paired-end experiments.")

parser.add_argument("--paired_std", dest="paired_stds", action='append',
                  help="Standard deviation for paired-end experiments.")



parser.add_argument("-o", "--output", dest="output_file", default="records.tsv",
                  help="The name of the output file.", metavar="file_name")

args = parser.parse_args()



def generate_file_md5(filename, blocksize=2**20):
    m = hashlib.md5()
    with open(filename, "rb" ) as f:
        while True:
            buf = f.read(blocksize)
            if not buf:
                break
            m.update( buf )
    return m.hexdigest()

if __name__ == '__main__':
    with open(args.output_file, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter=',')
        writer.writerow(["SERIES"])
        writer.writerow(["title", args.title])
        writer.writerow(["summary", args.summary])
        writer.writerow(["overall design", args.overall_design])
        if args.contributors is not None:
            for contributor in args.contributors:
                writer.writerow(["contributor", contributor])
        writer.writerow(["supplementary file", args.supplementary_file])
        writer.writerow(["SRA_center_name_code", args.sra_center_name_code])
        writer.writerow([])
        #The following section will be handeled dynamically since the headers also depend on the user input
        writer.writerow(["SAMPLES"])
        header_row = []
        header_row.append("Sample name")
        header_row.append("title")
        header_row.append("source name")
        header_row.append("organism")
        header_row.append("characteristics: cell type")
        header_row.append("characteristics: strain")
        header_row.append("characteristics: ChIP antibody")
        header_row.append("molecule")
        header_row.append("description")
        if args.samples_processed_input_files is not None:
            #need to get the max number of processed files for given sample here to be n
            #repeat the following n times
            max_count = 0;
            for i in range(len(args.samples_name)):
                curr_count = len(args.samples_processed_input_files[i].split(";"))
                if curr_count > max_count:
                    max_count = curr_count
            for x in range(max_count):
                header_row.append("processed data file")
        if args.samples_raw_input_files is not None:
            #need to get the max number of raw files for given sample here to be n
            #repeat the following n times
            max_count = 0;
            for i in range(len(args.samples_name)):
                curr_count = len(args.samples_raw_input_files[i].split(";"))
                if curr_count > max_count:
                    max_count = curr_count
            for x in range(max_count):
                header_row.append("raw file")
        writer.writerow(header_row)
        if args.samples_name is not None:
            row = []
            for i in range(len(args.samples_name)):
                row.append(args.samples_name[i])
                row.append(args.samples_title[i])
                row.append(args.samples_source_name[i])
                row.append(args.samples_organism[i])
                row.append(args.samples_cell_type[i])
                row.append(args.samples_strain[i])
                row.append(args.samples_chip[i])
                row.append(args.samples_molecule[i])
                row.append(args.samples_description[i])
                if args.samples_processed_input_files is not None:
                    my_list = args.samples_processed_input_files[i].split(";")
                    for item in my_list:
                        row.append(item)
                if args.samples_raw_input_files is not None:
                    my_list = args.samples_raw_input_files[i].split(";")
                    for item in my_list:
                        row.append(item)
                writer.writerow(row)
                row = []
        writer.writerow([])
        writer.writerow(["PROTOCOLS"])
        writer.writerow(["growth protocol", args.growth_protocol])
        writer.writerow(["treatment protocol", args.treatment_protocol])
        writer.writerow(["extract protocol", args.extract_protocol])
        writer.writerow(["library construction protocol", args.library_construction_protocol])
        writer.writerow(["library strategy", args.library_strategy])
        writer.writerow([])
        writer.writerow(["DATA PROCESSING PIPELINE"])
        if args.data_processing_steps is not None:
            for data_processing_step in args.data_processing_steps:
                writer.writerow(["data processing step", data_processing_step])
        writer.writerow(["genome build", args.genome_build])
        writer.writerow(["processed data files format and content", args.processed_data_files_format_content])
        writer.writerow([])
        writer.writerow(["PROCESSED DATA FILES"])
        writer.writerow(["file name", "file type", "file checksum"])
        if args.processed_files is not None:
            row = []
            for i in range(len(args.processed_files)):
                row.append(args.processed_file_names[i])
                row.append(os.path.splitext(args.processed_file_names[i])[1][1:])
                row.append(generate_file_md5(args.processed_files[i]))
                writer.writerow(row)
                row = []
        writer.writerow([])
        writer.writerow(["RAW FILES"])
        writer.writerow(["file name", "file type", "file checksum", "instrument model", "read length", "single or paired-end"])
        if args.raw_files is not None:
            row = []
            for i in range(len(args.raw_files)):
                row.append(args.raw_file_names[i])
                row.append(os.path.splitext(args.raw_file_names[i])[1][1:])
                row.append(generate_file_md5(args.raw_files[i]))
                row.append(args.raw_file_instrument_models[i])
                row.append(args.raw_file_read_lengths[i])
                row.append(args.raw_file_single_paired_choices[i])
                writer.writerow(row)
                row = []
        writer.writerow([])
        writer.writerow([])
        writer.writerow([])
        writer.writerow(["PAIRED-END EXPERIMENTS"])
        

        solid_experiments = False
        if args.paired_file_names1 is not None:
            row = []
            for i in range(len(args.paired_file_names1)):
                if args.paired_file_names3[i] != "None":
                    solid_experiments = True

        if (solid_experiments):
            writer.writerow(["file name 1", "file name 2", "file name 3", "file name 4", "average insert size", "standard deviation"])
        else:
            writer.writerow(["file name 1", "file name 2", "average insert size", "standard deviation"])

        if args.paired_file_names1 is not None:
            row = []
            for i in range(len(args.paired_file_names1)):
                row.append(args.paired_file_names1[i])
                row.append(args.paired_file_names2[i])
                if (solid_experiments):
                    row.append(args.paired_file_names3[i])
                    row.append(args.paired_file_names4[i])
                row.append(args.paired_avgs[i])
                row.append(args.paired_stds[i])
                writer.writerow(row)
                row = []

    tsvfile.close()