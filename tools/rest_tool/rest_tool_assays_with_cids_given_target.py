#!/usr/bin/env python


import argparse
import os
import sys

import readfile


def getIDofLine(line):
    arr = line.split(">")
    if len(arr) > 1:
        aid = arr[1].split("<")[0]
        return aid
    else:
        return "-1"


# get xml of all aids with cids for an activity
def getAllCidsForAssayActivity(activity):
    url = (
        "http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/activity/"
        + activity
        + "/aids/txt?list_return=listkey"
    )
    listkey = readfile.getresult(url)
    #    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/listkey/"+listkey+"/cids/xml"
    url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/25425,12345/cids/xml"
    print(("url: " + url))
    xml = readfile.getresult(url)

    # init parser
    handler = DictHandler()
    parser = sax.make_parser()
    parser.setContentHandler(handler)

    tempfile = open("tempfile", "w")
    # handle the last line, there is sometimes some random output
    lastline_arr = xml.split("\n")
    # print(lastline_arr)

    print ("l: ")
    print((len(lastline_arr)))
    lastline = lastline_arr[len(lastline_arr) - 1]
    print(("lastline: " + lastline))
    print(("lastline-2: " + lastline_arr[len(lastline_arr) - 2]))
    cidlastline = getIDofLine(lastline)
    aidkey = "-1"
    if cidlastline != "-1":
        i = len(lastline_arr) - 2
        # search for nex aid entry
        while i >= 0 and "AID" not in lastline_arr[i]:
            i -= 1
        if i >= 0:
            aid = getIDofLine(lastline_arr[i])
            if aid != "-1":
                aidkey = aid
    # remove the last line and put the array back together

    lastline_arr_list = list(lastline_arr)
    # lastline_arr_list.remove(lastline)
    xml2 = "\n".join(lastline_arr_list)
    tempfile.write(xml2)
    # add the last tags
    # tempfile.write("</Information></InformationList>")
    tempfile.close()
    parser.parse(open("tempfile", "r"))
    dic = handler.ergebnis

    # add the last line
    # if cidlastline != "-1":
    #    dic[aidkey].append(cidlastline)
    return dic


def main(args):
    aid_cid_dict = getAllCidsForAssayActivity(args.target)
    write_to_csv(aid_cid_dict, args.outfile)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--outfile", type=argparse.FileType("w"), help="Specify output file"
    )
    parser.add_argument("--target", type=str, help="Specify output file")
    if len(sys.argv) < 2:
        print("Too few arguments...")
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main(args)
