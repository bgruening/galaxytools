#!/usr/bin/env python
# Aufruf convert_graph.py --type type --operation op --id id --outformat format --outfile outfile

import sys, os
import networkx as nx
import argparse
import urllib2

#supported graph_types
output_types = ["tsv", "csv", "png", "json", "txt", "xml", "sdf", "asnt", "asnb", "jsonp"]


def getresult(url):
    try:
        connection = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return ""
    else:
        return connection.read().rstrip()
def main(args):
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/"+args.type+"/"
    if args.type == "assay":
        url+="aid/"
    elif args.type == "compound":
        url+="cid/"
    url+=args.id+"/"+args.operation+"/"+args.outformat
    print(url)
    print(args.type)
    data=getresult(url)
    file=args.outfile
    file.write(data)
    file.close()
    
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
    parser.add_argument('--outformat', type=str,
        help="Specify the format of the output", choices = output_types)
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        help="Specify one output file")
    if len(sys.argv) < 8:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
