#!/usr/bin/env python
"""
That tools tries so solve a simply mapping problem in the easiest and propably dirtiest way ever.
In galaxy we need often so called mapping of identifiers, e.g. EC-2-GO, UniProt-2-PDB, UniProt-2-PDB and so on.
Normally such mapping files are available in a tab-separated format.

AIM: Write a mapping tool, suitable for many different mapping-files with as 
less as possible modifications to the actual wrapper and with less admin-work as possible. 
So no database involved, no index creation ... just copy in the new mapping-file and it should work.

Idea: Writing a small wrapper around grep, take the information about the columns from a config-file stored next to the mapping file.
Inside the galaxy tool you ideally only need to specify the path to the mapping file.
Column headers and separaters can be stored in the config file.

"""

import sys, os
import argparse
import subprocess
import ConfigParser
import re, tempfile
import tempfile
from multiprocessing import Pool
import subprocess


def read_config_file(config_path="data.cfg"):
    with open(config_path) as configfile:
        for line in configfile:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            if line.lower().startswith("header="):
                header = line.split("=")[-1]
            if line.lower().startswith("separators="):
                sep = line.split("=")[-1].split()

    return (header, sep)


def chunks(list, size):
    """
    Yield successive n-sized chunks from list.
    """
    for i in xrange(0, len(list), size):
        yield list[i : i + size]


def single_grep(args):
    (pattern, tfile, mapping) = args
    subprocess.call(
        ["grep", "-w", "-i", "%s" % pattern, mapping], stdout=open(tfile, "w")
    )


def main(mapping, query, column_number, outfile):
    temp = tempfile.NamedTemporaryFile()

    (header, sep) = read_config_file(os.path.splitext(mapping)[0] + ".cfg")
    temp.write("%s\n" % header.replace("\\t", "\t"))
    temp.flush()

    iids = list()
    with open(query) as queryfile:
        for line in queryfile:
            line = line.strip()
            tokens = line.split("\t")
            iid = tokens[int(column_number) - 1]
            if iid.strip():
                # if the query has a DOT replace it with \.
                iids.append(re.escape(iid))

    patterns = list()
    temp_files = list()

    for chunck in chunks(iids, 10000):
        pattern = "\|".join(chunck)
        t = tempfile.NamedTemporaryFile(delete=False)
        t.close()
        patterns.append((pattern, t.name, mapping))
        temp_files.append(t.name)

    # single_grep(patterns[0])
    pool = Pool(processes=4)
    result = pool.map_async(single_grep, patterns)
    result.get()
    cmd = "cat %s > %s" % (" ".join(temp_files), outfile)
    os.system(cmd)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Small wrapper around grep to filter same mapping files."
    )

    parser.add_argument(
        "-m",
        "--mapping",
        dest="mapping_path",
        required=True,
        help="Path to the Mapping file.",
    )

    parser.add_argument(
        "-q",
        "--query",
        dest="query_path",
        required=True,
        help="Path to the query file that contains the ids.",
    )

    parser.add_argument(
        "-o",
        "--output",
        dest="output_path",
        required=True,
        help="Path to the output file.",
    )

    parser.add_argument(
        "-c",
        "--column-with-id",
        dest="column",
        default="c1",
        help="The column in which the id is stored.",
    )
    options = parser.parse_args()
    main(options.mapping_path, options.query_path, options.column, options.output_path)
