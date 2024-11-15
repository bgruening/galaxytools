__author__ = "Yasunobu OKAMURA"
__copyright__ = "Copyright (c) 2012 Y.Okamura"
__license__ = "GPL v3+"

import xml.parsers.expat

import networkx as nx


class XGMMLParserHelper(object):
    """"""

    def __init__(self, graph=nx.DiGraph()):
        """

        Arguments:
        - `graph`: Network X graph object
        """
        self._graph = graph
        self._parser = xml.parsers.expat.ParserCreate()
        self._parser.StartElementHandler = self._start_element
        self._parser.EndElementHandler = self._end_element
        self._tagstack = list()

        self._current_attr = dict()
        self._current_obj = dict()

    def _start_element(self, tag, attr):
        """

        Arguments:
        - `self`:
        - `tag`:
        - `attr`:
        """

        self._tagstack.append(tag)

        if tag == "node" or tag == "edge":
            self._current_obj = dict(attr)

        if tag == "att" and (
            self._tagstack[-2] == "node" or self._tagstack[-2] == "edge"
        ):
            if attr["type"] == "string":
                self._current_attr[attr["name"]] = attr["value"]
            elif attr["type"] == "real":
                self._current_attr[attr["name"]] = float(attr["value"])
            elif attr["type"] == "integer":
                self._current_attr[attr["name"]] = int(attr["value"])
            elif attr["type"] == "boolean":
                self._current_attr[attr["name"]] = bool(attr["value"])
            else:
                raise NotImplementedError(attr["type"])

    def _end_element(self, tag):
        """

        Arguments:
        - `self`:
        - `tag`:
        """

        if tag == "node":
            self._graph.add_node(
                self._current_obj["id"],
                label=self._current_obj["label"],
                **self._current_attr
            )
            # print 'add node', self._current_obj
        elif tag == "edge":
            self._graph.add_edge(
                self._current_obj["source"],
                self._current_obj["target"],
                **self._current_attr
            )

        self._tagstack.pop()

    def parseFile(self, file):
        """

        Arguments:
        - `self`:
        - `file`:
        """

        self._parser.ParseFile(file)

    def graph(self):
        """

        Arguments:
        - `self`:
        """

        return self._graph


def XGMMLWriter(file, graph, graph_name):
    """

    Arguments:
    - `graph`:
    """

    print("""<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<graph directed="1"  xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#" xmlns="http://www.cs.rpi.edu/XGMML">
<att name="selected" value="1" type="boolean" />
<att name="name" value="{0}" type="string"/>
<att name="shared name" value="{0}" type="string"/>
""".format(
        graph_name
    ), file=file)

    for onenode in graph.nodes(data=True):
        id = onenode[0]
        attr = dict(onenode[1])

        if "label" in attr:
            label = attr["label"]
            del attr["label"]
        else:
            label = id

        print('<node id="{id}" label="{label}">'.format(id=id, label=label), file=file)
        for k, v in attr.iteritems():
            print('<att name="{}" value="{}" type="string" />'.format(k, v), file=file)
        print("</node>", file=file)

    for oneedge in graph.edges(data=True):
        print('<edge source="{}" target="{}">'.format(oneedge[0], oneedge[1]), file=file)
        for k, v in oneedge[2].iteritems():
            print('<att name="{}" value="{}" type="string" />'.format(k, v), file=file)
        print("</edge>", file=file)
    print("</graph>", file=file)
