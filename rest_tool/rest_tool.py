#!/usr/bin/env python
# Aufruf convert_graph.py --type type --operation op --id id --outfile outfile

import sys, os
import argparse

import readfile

txt_output=["cids", "summary", "synonyms" ]

def main(args):
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/"+args.type+"/"
    if args.type == "assay":
        url+="aid/"
    elif args.type == "compound":
        url+="cid/"
    elif args.type == "substance":
        url+="sid/"
    if args.idfile is None:
        idstring=str(args.id)
    else:
        idlist=readfile.getListFromFile(args.idfile)
        idstring=",".join(idlist)
    url+=idstring+"/"+args.operation+"/"
    if args.operation in txt_output:
        url+="txt"
    else:
        url+="csv"
    print(url)
    data=readfile.getresult(url)
    outfile=args.outfile
    outfile.write(data)
    outfile.close()
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--type', type=str,
        help="That you want BioAssay Compund ...")
    parser.add_argument('--id', type=str,
        help="Specify the ID")
    parser.add_argument('--operation', type=str,
        help="Specify the operation")
    parser.add_argument('--property-value', type=str,
        help="Specify the property")
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        help="Specify one output file")
    parser.add_argument('--idfile', type=argparse.FileType('r'),
        help="Specify a file with a list of ids, one per line")
    if len(sys.argv) < 8:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
