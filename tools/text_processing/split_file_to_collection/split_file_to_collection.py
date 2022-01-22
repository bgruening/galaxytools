#!/usr/bin/env python

import argparse
import math
import os
import random
import re

# configuration of the splitting for specific file types
# - regular expression matching the record separator ('' if not splitting by regex but by number of lines)
# - number of lines to split after (0 if not splitting by number of lines but regex)
# - a boolean indicating if the record separator is at the end of the record
#
# new file types can be added by appending to this dict,
# updating the parser, and adding a new type option in the Galaxy wrapper
FILETYPES = {
    "fasta": (r"^>", 0, False),
    "fastq": (r"", 4, False),
    "tabular": (r"", 1, False),
    "txt": (r"", 1, False),
    "mgf": (r"^BEGIN IONS", 0, False),
    "sdf": (r"\$\$\$\$", 0, True),
}


def main():
    ps = parser_cli()
    args = vars(ps.parse_args())

    # get args and validate
    in_file = args["in"]
    if not os.path.isfile(args["in"]):
        raise FileNotFoundError("Input file does not exist")

    out_dir = args["out_dir"]
    if not os.path.isdir(args["out_dir"]):
        raise FileNotFoundError("out_dir is not a directory")

    top = args["top"]
    if top < 0:
        raise ValueError("Number of header lines cannot be negative")

    ftype = args["ftype"]

    assert (
        ftype != "generic" or args["generic_re"] is not None
    ), "--generic_re needs to be given for generic input"

    if args["ftype"] == "tabular" and args["by"] == "col":
        args["match"] = replace_mapped_chars(args["match"])
        args["sub"] = replace_mapped_chars(args["sub"])
        split_by_column(args, in_file, out_dir, top)
    else:
        args["generic_re"] = replace_mapped_chars(args["generic_re"])
        split_by_record(args, in_file, out_dir, top, ftype)


def parser_cli():
    parser = argparse.ArgumentParser(
        description="split a file into multiple files. "
        + "Can split on the column of a tabular file, "
        + "with custom and useful names based on column value."
    )
    parser.add_argument("--in", "-i", required=True, help="The input file")
    parser.add_argument(
        "--out_dir",
        "-o",
        default=os.getcwd(),
        help="The output directory",
        required=True,
    )
    parser.add_argument(
        "--file_names",
        "-a",
        help="If not splitting by column, the base name of the new files",
    )
    parser.add_argument(
        "--file_ext",
        "-e",
        help="If not splitting by column,"
        + " the extension of the new files (without a period)",
    )
    parser.add_argument(
        "--ftype",
        "-f",
        help="The type of the file to split",
        required=True,
        choices=["mgf", "fastq", "fasta", "sdf", "tabular", "txt", "generic"],
    )
    parser.add_argument(
        "--by",
        "-b",
        help="Split by line or by column (tabular only)",
        default="row",
        choices=["col", "row"],
    )
    parser.add_argument(
        "--top",
        "-t",
        type=int,
        default=0,
        help="Number of header lines to carry over to new files.",
    )
    parser.add_argument(
        "--rand",
        "-r",
        help="Divide records randomly into new files",
        action="store_true",
    )
    parser.add_argument(
        "--seed",
        "-x",
        help="Provide a seed for the random number generator. "
        + 'If not provided and args["rand"]==True, then date is used',
        type=int,
    )
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--numnew",
        "-n",
        type=int,
        default=1,
        help="Number of output files desired. Not valid for splitting on a column. Not compatible with chunksize and will be ignored if both are set.",
    )
    group.add_argument(
        "--chunksize",
        "-k",
        type=int,
        default=0,
        help="Number of records by file. Not valid for splitting on a column",
    )
    parser.add_argument(
        "--batch",
        action="store_true",
        help="Distribute files to collection while maintaining order. Ignored if splitting on column.",
    )
    generic = parser.add_argument_group("Arguments controling generic splitting")
    group = generic.add_mutually_exclusive_group()
    group.add_argument(
        "--generic_re",
        "-g",
        default="",
        help="Regular expression indicating the start of a new record (only for generic)",
        required=False,
    )
    group.add_argument(
        "--generic_num",
        type=int,
        default=0,
        help="Length of records in number of lines (only for generic)",
        required=False,
    )
    generic.add_argument(
        "--split_after",
        "-p",
        action="store_true",
        help="Split between records after separator (default is before). "
        + "Only for generic splitting by regex - specific ftypes are always split in the default way",
    )
    bycol = parser.add_argument_group("If splitting on a column")
    bycol.add_argument(
        "--match",
        "-m",
        default="(.*)",
        help="The regular expression to match id column entries",
    )
    bycol.add_argument(
        "--sub",
        "-s",
        default=r"\1",
        help="The regular expression to substitute in for the matched pattern.",
    )
    bycol.add_argument(
        "--id_column",
        "-c",
        default="1",
        help="Column that is used to name output files. Indexed starting from 1.",
        type=int,
    )
    return parser


def replace_mapped_chars(pattern):
    """
    handles special escaped characters when coming from galaxy
    """
    mapped_chars = {"'": "__sq__", "\\": "__backslash__"}
    for key, value in mapped_chars.items():
        pattern = pattern.replace(value, key)
    return pattern


def split_by_record(args, in_file, out_dir, top, ftype):
    # get configuration (record separator, start at end) for given filetype
    sep, num, sep_at_end = FILETYPES.get(
        ftype, (args["generic_re"], args["generic_num"], args["split_after"])
    )
    sep = re.compile(sep)

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

    # determine
    # - the number of records that should be stored per file
    #   (done always, even if used only for batch mode)
    # - if the separator is a the start / end of the record
    n_per_file = math.inf
    if (
        chunksize != 0 or batch
    ):  # needs to be calculated if either batch or chunksize are selected
        with open(in_file) as f:
            # read header lines
            for i in range(top):
                f.readline()
            n_records = 0
            for line in f:
                if (num == 0 and re.match(sep, line) is not None) or (
                    num > 0 and n_records % num == 0
                ):
                    n_records += 1
                    last_line_matched = True
                else:
                    last_line_matched = False
            if sep_at_end and not last_line_matched:
                n_records += 1

        # if there are fewer records than desired files
        numnew = min(numnew, n_records)
        # approx. number of records per file
        if chunksize == 0:  # i.e. no chunking
            n_per_file = n_records // numnew
        else:
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

    newfile_names = [
        os.path.join(out_dir, "%s_%06d%s" % (new_file_base[0], count, new_file_base[1]))
        for count in range(0, numnew)
    ]
    # bunch o' counters
    # index to list of new files
    if rand:
        new_file_counter = int(math.floor(random.random() * numnew))
    else:
        new_file_counter = 0
    new_file = open(newfile_names[new_file_counter], "a")
    # to contain header specified by top
    header = ""
    # keep track of the files that have been opened so far
    fresh_files = set(range(numnew))

    # keep track in loop of number of records in each file
    # only used in batch
    records_in_file = 0

    # open file
    with open(in_file, "r") as f:
        # read header
        for i in range(top):
            header += f.readline()

        record = ""
        for line_no, line in enumerate(f):
            # check if beginning of line is record sep
            # if beginning of line is record sep, either start record or finish one
            if (num == 0 and re.match(sep, line) is not None) or (
                num > 0 and line_no % num == 0
            ):
                # this only happens first time through
                if record == "":
                    record += line
                else:
                    # if is in fresh_files, write header and drop from freshFiles
                    if new_file_counter in fresh_files:
                        new_file.write(header)
                        fresh_files.remove(new_file_counter)

                    if sep_at_end:
                        record += line
                    # write record to file
                    new_file.write(record)
                    if not sep_at_end:
                        record = line
                    else:
                        record = ""

                    # change destination file
                    if rand:
                        new_file_counter = int(math.floor(random.random() * numnew))
                        new_file.close()
                        new_file = open(newfile_names[new_file_counter], "a")
                    elif batch:
                        # number of records read per file
                        records_in_file += 1
                        # have we reached the max for each file?
                        # if so, switch file
                        if records_in_file >= n_per_file:
                            new_file_counter = (new_file_counter + 1) % numnew
                            records_in_file = 0  # reset to 0
                            new_file.close()
                            new_file = open(newfile_names[new_file_counter], "a")
                    else:
                        new_file_counter = (new_file_counter + 1) % numnew
                        new_file.close()
                        new_file = open(newfile_names[new_file_counter], "a")
            # if beginning of line is not record sep, we must be inside a record
            # so just append
            else:
                record += line
        # after loop, write final record to file
        new_file.write(record)
        new_file.close()


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
    files = set()

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
            fields = re.split(r"\t", line.strip("\n"))

            # get id column value
            id_col_val = fields[id_col]

            # use regex to get new file name
            out_file_name = re.sub(match, sub, id_col_val)
            out_file_path = os.path.join(out_dir, out_file_name)

            # write
            with open(out_file_path, "a") as current_new_file:
                if out_file_name not in files:
                    current_new_file.write(header)
                    files.add(out_file_name)
                current_new_file.write(line)


if __name__ == "__main__":
    main()
