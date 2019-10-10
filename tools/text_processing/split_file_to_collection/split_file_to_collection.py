#!/usr/bin/env python

import argparse
import os
import re
import random
import math


"""
regexes that indicate the *beginning* of a record
new file types can be added by appending to this dict,
updating the parser, and adding a new type option in the Galaxy wrapper
"""
FILETYPES = {'fasta': '^>',
             'fastq': '^@',
             'tabular': '^.*',
             'txt': '^.*',
             'mgf': '^BEGIN IONS',
             'sdf': '\$\$\$\$',
             }


def main():
    ps = parser_cli()
    args = vars(ps.parse_args())

    # get args and validate
    in_file = args["in"]
    if not os.path.isfile(args["in"]):
        raise FileNotFoundError('Input file does not exist')

    out_dir = args["out_dir"]
    if not os.path.isdir(args["out_dir"]):
        raise FileNotFoundError('out_dir is not a directory')

    top = args["top"]
    if top < 0:
        raise ValueError("Number of header lines cannot be negative")

    ftype = args["ftype"]

    assert ftype != "generic" or args["generic_re"] != None, "--generic_re needs to be given for generic input"

    if args["ftype"] == "tabular" and args["by"] == "col":
        args["match"] = replace_mapped_chars(args["match"])
        args["sub"] = replace_mapped_chars(args["sub"])
        split_by_column(args, in_file, out_dir, top)

    else:
        split_by_record(args, in_file, out_dir, top, ftype)


def parser_cli():
    parser = argparse.ArgumentParser(description="split a file into multiple files. " +
                                                 "Can split on the column of a tabular file, " +
                                                 "with custom and useful names based on column value.")
    parser.add_argument('--in', '-i', required=True, help="The input file")
    parser.add_argument('--out_dir', '-o', default=os.getcwd(), help="The output directory", required=True)
    parser.add_argument('--file_names', '-a', help="If not splitting by column, the base name of the new files")
    parser.add_argument('--file_ext', '-e', help="If not splitting by column," +
                                                 " the extension of the new files (without a period)")
    parser.add_argument('--ftype', '-f', help="The type of the file to split", required = True,
        choices=["mgf", "fastq", "fasta", "sdf", "tabular", "txt", "generic"])
    parser.add_argument('--generic_re', '-g', help="Regular expression indicating the start of a new record (only for generic)", required = False)
    parser.add_argument('--by', '-b', help="Split by line or by column (tabular only)",
        default = "row", choices = ["col", "row"])
    parser.add_argument('--top', '-t', type=int, default=0, help="Number of header lines to carry over to new files. " +
                                                                 "(tabular only).")
    parser.add_argument('--rand', '-r', help="Divide records randomly into new files", action='store_true')
    parser.add_argument('--seed', '-x', help="Provide a seed for the random number generator. " +
                                             "If not provided and args[\"rand\"]==True, then date is used", type=int)
    parser.add_argument('--numnew', '-n', type=int, default = 1,
                        help="Number of output files desired. Not valid for splitting on a column. Not compatible with chunksize and will be ignored if both are set.")
    parser.add_argument('--chunksize', '-k', type=int, default = 0,
                        help="Number of records by file. Not valid for splitting on a column")
    parser.add_argument('--batch', action='store_true',
                        help="Distribute files to collection while maintaining order. Ignored if splitting on column.")
    parser.add_argument('--split_after', '-p', action='store_true',
                        help="Split between records after separator (default is before)." + 
                         "Only for generic - specific ftypes are always split in the default way")
    bycol = parser.add_argument_group('If splitting on a column')
    bycol.add_argument('--match', '-m', default = "(.*)", help="The regular expression to match id column entries")
    bycol.add_argument('--sub', '-s', default = r'\1',
                       help="The regular expression to substitute in for the matched pattern.")
    bycol.add_argument('--id_column', '-c', default="1",
                       help="Column that is used to name output files. Indexed starting from 1.", type=int)
    return parser


def close_files(file_list):
    # finally, close all files
    for open_file in file_list:
        open_file.close()


def replace_mapped_chars(pattern):
    """
    handles special escaped characters when coming from galaxy
    """
    mapped_chars = {'\'': '__sq__', '\\': '__backslash__'}
    for key, value in mapped_chars.items():
        pattern = pattern.replace(value, key)
    return pattern


def split_by_record(args, in_file, out_dir, top, ftype):
    # get record separator for given filetype
    sep = re.compile(FILETYPES.get(ftype, args["generic_re"]))

    chunksize = args["chunksize"]
    numnew = args["numnew"]

    # random division
    rand = args["rand"]
    seed = args["seed"]
    if seed:
        random.seed(seed)
    else:
        random.seed()

    # batched division (maintains order)
    batch = args["batch"]

    
    if chunksize != 0 or batch: # needs to be calculated if either batch or chunksize are selected
        # define n_per_file so we don't get a warning about ref before assignment
        n_per_file = math.inf

        # number of records
        with open(in_file) as f:
            i = 0
            for line in f:
                if re.match(sep, line) is not None:
                    i+=1
            n_records = i + 1
        if top:
            n_records -= top  # don't count the top lines
        
        if chunksize == 0: # i.e. no chunking
            # approx. number of lines per file
            n_per_file = n_records // numnew
        else:
            # approx. number of lines per file
            numnew = n_records // chunksize
            n_per_file = chunksize




    # make new files
    # strip extension of old file and add number
    custom_new_file_name = args["file_names"]
    custom_new_file_ext = "." + args["file_ext"]
    if custom_new_file_name is None:
        new_file_base = os.path.splitext(os.path.basename(in_file))
    else:
        new_file_base = [custom_new_file_name, custom_new_file_ext]

    newfiles = [
        open(os.path.join(out_dir, "%s_%06d%s" % (new_file_base[0], count, new_file_base[1])) , "w")
        for count in range(0, numnew)
    ]

    # bunch o' counters
    # index to list of new files
    new_file_counter = 0

    # used for top
    # number of lines read so far
    n_read = 0
    # to contain header specified by top
    header = ""
    # keep track of the files that have been opened so far
    fresh_files = {i for i in range(0, numnew)}

    # keep track in loop of number of records in each file
    # only used in batch
    records_in_file = 0

    # open file
    with open(in_file, "r") as file:
        record = ""
        for line in file:
            n_read += 1
            if n_read <= top:
                header += line
                continue
            # check if beginning of line is record sep
            # if beginning of line is record sep, either start record or finish one
            if re.match(sep, line) is not None:
                # this only happens first time through
                if record == "":
                    record += line
                else:
                    # if is in fresh_files, write header and drop from freshFiles
                    if new_file_counter in fresh_files:
                        newfiles[new_file_counter].write(header)
                        fresh_files.remove(new_file_counter)
                    
                    if ftype != "sdf" and args["split_after"] == False:
                        # write record to file
                        newfiles[new_file_counter].write(record)

                        # if not the first time through, we assign the new record
                        record = line
                                                
                    else:  # for sdf we want to write the line to the record before starting a new one
                        record += line
                        newfiles[new_file_counter].write(record)
                        record = ""
                        
                    # change destination file
                    if rand:
                        new_file_counter = int(math.floor(random.random() * numnew))
                    elif batch:
                        # number of records read per file
                        records_in_file += 1
                        # have we reached the max for each file?
                        # if so, switch file
                        if records_in_file >= n_per_file:
                            new_file_counter = (new_file_counter + 1) % numnew
                            records_in_file = 0  # reset to 0
                    else:
                        new_file_counter = (new_file_counter + 1) % numnew
            # if beginning of line is not record sep, we must be inside a record
            # so just append
            else:
                record += line
        # after loop, write final record to file
        newfiles[new_file_counter].write(record)
    # close new files
    close_files(newfiles)


def split_by_column(args, in_file, out_dir, top):

    # shift to 0-based indexing
    id_col = int(args["id_column"]) - 1

    try:
        match = re.compile(args["match"])
    except re.error:
        print("ERROR: Match (-m) supplied is not valid regex.")
        raise

    sub = args["sub"]

    # set of file names
    new_files = dict()

    # keep track of how many lines have been read
    n_read = 0
    header = ""
    with open(in_file) as file:
        for line in file:
            # if still in top, save to header
            n_read += 1
            if n_read <= top:
                header += line
                continue
            # split into columns, on tab
            fields = re.split(r'\t', line.strip('\n'))

            # get id column value
            id_col_val = fields[id_col]

            # use regex to get new file name
            out_file_name = re.sub(match, sub, id_col_val)
            out_file_path = os.path.join(out_dir, out_file_name)

            # write
            if out_file_name not in new_files.keys():
                # open file (new, so not already open)
                current_new_file = open(out_file_path, "w")
                current_new_file.write(header)
                current_new_file.write(line)
                # add to dict
                new_files[out_file_name] = current_new_file
            else:
                # file is already open, so just write to it
                new_files[out_file_name].write(line)

    # finally, close all files
    close_files(new_files.values())


if __name__ == "__main__":
    main()
