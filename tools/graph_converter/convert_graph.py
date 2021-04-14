#!/usr/bin/env python
# Aufruf convert_graph.py --infile datei --informat typ --outfile ausgabedatei --outformat ausgabetyp

import sys, os
import networkx as nx
import argparse
import json

from xgmml_networkx import XGMMLParserHelper, XGMMLWriter
from networkx.readwrite import json_graph

# supported graph_types
graph_types = ["gml", "yaml", "gspan", "xgmml", "gexf", "graphml", "json", "pajek"]

func_dic_read = {
    "gml": nx.read_gml,
    "yaml": nx.read_yaml,
    "gexf": nx.read_gexf,
    "graphml": nx.read_graphml,
    "pajek": nx.read_pajek,
}

func_dic_write = {
    "gml": nx.write_gml,
    "yaml": nx.write_yaml,
    "gexf": nx.write_gexf,
    "graphml": nx.write_graphml,
    "pajek": nx.write_pajek,
}

# completely supported types by networkx
completely_supported_types = ["gml", "gexf", "yaml", "graphml", "pajek"]


def read_gspan(infile):
    G = nx.DiGraph()
    idoffset = 0
    old_id_start = 0
    for line in infile:
        line_split = line.split(" ")
        length_split = len(line_split)
        if line[0] == "v":
            G.add_node(idoffset, label=line_split[2].strip())
            idoffset += 1
        elif line[0] == "e":
            if length_split < 3:
                raise InvalidGraph(line)
            elif length_split > 3:
                G.add_edge(
                    old_id_start + int(line_split[1]),
                    old_id_start + int(line_split[2]),
                    label=line_split[3].strip(),
                )
            else:
                G.add_edge(
                    old_id_start + int(line_split[1]),
                    old_id_start + int(line_split[2]),
                    label="",
                )
        elif line[0] == "t":
            # its a new subgraph
            # idoffset*=1
            old_id_start = idoffset
    # print(nx.is_connected(G))
    return G


def write_gspan(graph, outfile):
    # get all subgraphs only works with undirected
    subgraphs = nx.connected_components(graph.to_undirected())
    id_count = 1
    node_count = 0
    # get labels
    label_dic = nx.get_node_attributes(graph, "label")
    for s in subgraphs:
        node_count_tree = 0
        node_dict = {}
        outfile.write("t # id " + str(id_count) + "\n")
        # for every node in subgraph
        for v in sorted(s):
            # node id restart from 0 for every sub graph
            node_dict[v] = node_count_tree
            outfile.write("v " + str(node_count_tree) + " " + label_dic[v] + " \n")
            node_count_tree += 1
            node_count += 1

        # all edges adjacent to a node of s
        edges = nx.edges(graph, s)
        for e in sorted(edges):
            # print(graph[e[0]][e[1]])
            try:
                outfile.write(
                    "e "
                    + str(node_dict[e[0]])
                    + " "
                    + str(node_dict[e[1]])
                    + " "
                    + graph[e[0]][e[1]]["label"]
                    + "\n"
                )
            except KeyError:
                outfile.write("e " + str(node_dict[e[0]]) + " " + str(node_dict[e[1]]))

        id_count += 1


def read_json(file):
    json_string = file.read()
    # print(json_string)
    json_dict = json_graph.loads(json_string)
    # print(json_dict)
    # return json_graph.node_link_graph(json_dict, True, False)
    return json_dict


def write_json(graph, outfile):
    json_dict = json_graph.node_link_data(graph)
    json_string = json_graph.dumps(json_dict)
    outfile.write(json_string)
    # print("did it")


def main(args):

    if args.informat not in graph_types:
        print "EXCEPTION COMPUTER EXPLODING"
    # everything networkx can do by itself ;)
    elif args.informat in completely_supported_types:
        function = func_dic_read[args.informat]
        graph = function(args.infile)
    elif args.informat == "gspan":
        graph = read_gspan(args.infile)
    elif args.informat == "json":
        graph = read_json(args.infile)
    elif args.informat == "xgmml":
        xgmml = XGMMLParserHelper()
        xgmml.parseFile(args.infile)
        graph = xgmml.graph()

    if args.outformat in completely_supported_types:
        function = func_dic_write[args.outformat]
        function(graph, args.outfile)
    elif args.outformat == "gspan":
        write_gspan(graph, args.outfile)
    elif args.outformat == "json":
        write_json(graph, args.outfile)
    elif args.outformat == "xgmml":
        # xgmml=XGMMLParserHelper(graph)
        # xgmml.parseFile(open(sys.argv[1]))
        a = XGMMLWriter(args.outfile, graph, "MyGraph")


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--infile",
        type=argparse.FileType("r"),
        help="Specify the input file representing a graph",
    )
    parser.add_argument(
        "--outfile", type=argparse.FileType("w"), help="Specify one output file"
    )
    parser.add_argument(
        "--informat",
        type=str,
        help="Specify the format of the input graph",
        choices=graph_types,
    )
    parser.add_argument(
        "--outformat",
        type=str,
        help="Specify the format of the output graph",
        choices=graph_types,
    )
    if len(sys.argv) < 8:
        print "Too few arguments..."
        parser.print_help()
        exit(1)
    args = parser.parse_args()
    main(args)
