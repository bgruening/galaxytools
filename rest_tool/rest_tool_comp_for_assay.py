#!/usr/bin/env python
# Aufruf convert_graph.py --aid list of ids --aid-from-file file

import sys, os
import networkx as nx
import argparse
import urllib2, urllib, httplib

import readfile
#supported graph_types
#output_types = ["tsv", "csv", "png", "json", "txt", "xml", "sdf", "asnt", "asnb", "jsonp"]

        
#get the cids for bioassay aid
def getCompoundList(aidlist):
    aidliststring=",".join(aidlist)
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+aidliststring+"/cids/txt"
    data=readfile.getresult(url)
    return data
        
def main(args):
    if args.aidfile is None:
        aidlist=args.aid.split(",")
    else:
        aidlist=readfile.getListFromFile(args.aidfile)
    cids=getCompoundList(aidlist)
    args.outfile.write(cids)
    args.outfile.close()
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--aid', type=str,
        help="AIDs of the BioAssay")
    parser.add_argument('--aidfile', type=argparse.FileType('r'),
        help="Specify a file with a list of aids, one per line")
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        help="Specify output file")
    if len(sys.argv) < 2:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
