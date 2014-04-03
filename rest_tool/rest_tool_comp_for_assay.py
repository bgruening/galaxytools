#!/usr/bin/env python

import sys, os
import argparse
import tempfile
import readfile
import rest_tool_functions

        
#get the cids for bioassay aid
def get_aid_cid_dict_for_list(aidlist):
    aidliststring=",".join(aidlist)
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+aidliststring+"/cids/xml"
    xml=readfile.getresult(url)
    tmp = tempfile.TemporaryFile() 
    tmp.write(xml)
    tmp.seek(0)
    dic=rest_tool_functions.give_aid_cid_dict_from_xml(tmp)
    tmp.close()
    return dic
        
def main(args):
    if args.aid_file is None:
        aidlist=args.aid.split(",")
    else:
        aidlist=readfile.getListFromFile(args.aid_file)
    dic=get_aid_cid_dict_for_list(aidlist)
    rest_tool_functions.write_to_csv(dic, args.outfile)
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--aid', type=str,
        help="AIDs of the BioAssay")
    parser.add_argument('--aid-file', dest="aid_file", type=argparse.FileType('r'),
        help="Specify a file with a list of aids, one per line")
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        help="Specify output file")
    if len(sys.argv) < 2:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
