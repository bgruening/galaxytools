#!/usr/bin/env python
# Aufruf convert_graph.py --aid list of ids --aid-from-file file

import sys, os
import networkx as nx
import argparse
import urllib2, urllib, httplib

#supported graph_types
#output_types = ["tsv", "csv", "png", "json", "txt", "xml", "sdf", "asnt", "asnb", "jsonp"]


def getresult(url):
    try:
        connection = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return ""
    else:
        return connection.read().rstrip()
        
        
def getresultPost(url, aidList):
    values={'aids': ",".join(aidList)}
    data=urllib.urlencode(values)
    print(data)
    
    headers = {"Content-type": "application/x-www-form-urlencoded"}
    conn = httplib.HTTPConnection("pubchem.ncbi.nlm.nih.gov")
    conn.request("POST", "/rest/pug/assay/aid/aids/csv", data, headers)
    response = conn.getresponse()
    conn.close()
    return response.read()


        
#get the cids for bioassay aid
def getCompoundList(aid):
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/aid/"+aid+"/cids/txt"
    data=getresult(url)
    return data
        
def main(args):
    cids=getCompoundList(args.aid)
    args.outfile.write(cids)
    args.outfile.close()
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--aid', type=str,
        help="AIDs of the BioAssay")
    parser.add_argument('--aid-from-file', type=argparse.FileType('r'),
        help="Specify a file with a list of aids, one per line")
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        help="Specify output file")
    if len(sys.argv) < 2:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
