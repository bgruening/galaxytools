#!/usr/bin/env python

import argparse
import os
import sys

import readfile

txt_output = ["cids", "aids", "sids", "synonyms"]
csv_output = ["assaysummary", "property"]
check_for_id_type = ["cids", "aids", "sids"]


def main(args):
    url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/" + args.type + "/"
    if args.type == "assay":
        url += "aid/"
    elif args.type == "compound":
        url += "cid/"
    elif args.type == "substance":
        url += "sid/"
    if args.id_file is None:
        idstring = str(args.id)
    else:
        idlist = readfile.getListFromFile(args.id_file)
        idstring = ",".join(idlist)
    url += idstring + "/" + args.operation + "/"
    if args.operation == "property":
        url += args.property_value + "/"
    if args.operation in csv_output:
        url += "csv"
    elif args.operation in txt_output:
        url += "txt"
    else:
        url += "xml"
    if args.operation in check_for_id_type and not args.id_type is None:
        url += "?" + args.operation + "_type=" + args.id_type
    print(("The constructed REST URL is: %s" % url))
    data = readfile.getresult(url)
    outfile = args.outfile
    outfile.write(data)
    outfile.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--type", type=str, required=True, help="That you want BioAssay Compund ..."
    )
    parser.add_argument("--id", type=str, help="Specify the ID")
    parser.add_argument(
        "--operation", type=str, required=True, help="Specify the operation"
    )
    parser.add_argument(
        "--property-value", dest="property_value", type=str, help="Specify the property"
    )
    parser.add_argument(
        "--id-type", dest="id_type", type=str, help="Specify the property"
    )
    parser.add_argument(
        "--outfile",
        type=argparse.FileType("w"),
        required=True,
        help="Specify one output file",
    )
    parser.add_argument(
        "--id-file",
        dest="id_file",
        type=argparse.FileType("r"),
        help="Specify a file with a list of ids, one per line",
    )
    if len(sys.argv) < 8:
        print("Too few arguments...")
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main(args)
