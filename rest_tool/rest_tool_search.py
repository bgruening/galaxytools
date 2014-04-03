#!/usr/bin/env python

import sys, os
import argparse
import readfile

def main(args):
    #search for acitivity or target
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/"+args.type+"/name/"+args.name
    if args.type == "assay":
        url+="/aids"
    elif args.type == "compound":
        url+="/cids"
    else:
        url+="/sids"
    url+="/txt"
    #print("url: "+url)
    data=readfile.getresult(url)
    args.outfile.write(data)
    args.outfile.close()
    
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--name', type=str, required=True,
        help="Enter the name")
    parser.add_argument('--type', type=str, required=True,
        help="What you want to search for")
    parser.add_argument('--outfile', type=argparse.FileType('w'), required=True,
        help="Specify output file")
    if len(sys.argv) < 2:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
