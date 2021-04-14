#!/usr/bin/env python

import argparse

import readfile
import rest_tool_functions

"""
Get all DICs from belonging to bioassay IDs (AIDs)
"""


def main(args):
    if args.aid_file is None:
        aidlist = args.aid.split(",")
    else:
        aidlist = readfile.getListFromFile(args.aid_file)
    aidliststring = ",".join(aidlist)
    url = (
        "http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"
        + aidliststring
        + "/cids/xml"
    )
    print("The constructed REST URL is: %s" % url)
    dic = rest_tool_functions.get_dict_key_value(url, "AID", "CID")
    rest_tool_functions.write_to_sf(dic, args.outfile, "\t")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--aid", type=str, dest="aid", help="AIDs of the BioAssay")
    parser.add_argument(
        "--aid-file",
        dest="aid_file",
        type=argparse.FileType("r"),
        help="Specify a file with a list of aids, one per line",
    )
    parser.add_argument(
        "--outfile", type=argparse.FileType("w"), help="Specify output file"
    )

    args = parser.parse_args()
    main(args)
