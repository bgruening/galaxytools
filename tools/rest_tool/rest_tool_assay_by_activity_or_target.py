#!/usr/bin/env python

import argparse
import os
import sys

import readfile


def main(args):
    # search for acitivity or target
    url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/"
    if args.activity is None:
        # target
        url += "target/%s/%s" % (args.targettype, args.targetid)
    else:
        url += "activity/" + args.activity
    url += "/aids/txt"
    data = readfile.getresult(url)
    args.outfile.write(data)
    args.outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--activity", type=str, help="Activities you are looking for")
    parser.add_argument(
        "--target-type", dest="target_type", type=str, help="The target identifier type"
    )
    parser.add_argument(
        "--target-id", dest="target_id", type=str, help="The specific target"
    )
    parser.add_argument(
        "--outfile",
        type=argparse.FileType("w"),
        required=True,
        help="Specify output file",
    )
    if len(sys.argv) < 2:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main(args)
