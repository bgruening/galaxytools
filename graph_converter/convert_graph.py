#!/usr/bin/env python
# Aufruf convert_graph.py --infile datei --informat typ --outfile ausgabedatei --outformat ausgabetyp

import sys, os
import networkx as nx
import argparse
from xgmml_networkx import XGMMLParserHelper, XGMMLWriter


#supported graph_types
graph_types = ["gml", "yaml", "gspan", "xgmml", "gexf", "graphml", "json", "pajek"]

#completely supported types by networkx
complete_supported_types = ["gml", "gexf", "yaml", "graphml", "pajek"]

def main( args ):

    if args.informat not in graph_types:
        print "EXCEPTION COMPUTER EXPLODING"
    elif args.informat in complete_supported_types:
        function="nx.read_"+args.informat
        function(args.infile)
    """
    if args.informat == "gml":
        graph=nx.read_gml(args.infile)  
    elif args.informat == "gexf":
        graph=nx.read_gexf(args.infile)
    elif args.informat == "yaml":
        graph=nx.read_yaml(args.infile)
    elif args.informat == "xgmml":
        xgmml=XGMMLParserHelper()
        xgmml.parseFile(args.infile)
        graph=xgmml.graph()
        print(graph)
    # everything networkx can do by itself ;)
    """
    if args.outformat == "gml":
        nx.write_gml(graph, args.outfile)
    elif args.outformat == "gexf":
        nx.write_gexf(graph, args.outfile)
    elif args.outformat == "yaml":
        nx.write_yaml(graph, args.outfile)
    elif args.outformat == "xgmml":  
    #xgmml=XGMMLParserHelper(graph)
    #xgmml.parseFile(open(sys.argv[1]))
        a=XGMMLWriter(args.outfile, graph, "MyGraph")
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('--infile', type=argparse.FileType('r'),
        help="Specify the input file representing a graph")
    parser.add_argument('--outfile', type=argparse.FileType('w'),
        help="Specify one output file")
    parser.add_argument('--informat', type=str,
        help="Specify the format of the input graph", choices = graph_types)
    parser.add_argument('--outformat', type=str,
        help="Specify the format of the output graph", choices = graph_types)
    if len(sys.argv) < 4:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main( args )
