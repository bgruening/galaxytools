#!/usr/bin/env python
"""
uniprot python interface
to access the uniprot database

Based on work from Jan Rudolph: https://github.com/jdrudolph/uniprot
available services:
    map
    retrieve
"""

import requests
import sys, argparse

url = 'http://www.uniprot.org/'

def _retrieve(query, format='txt'):
    """_retrieve is not meant for use with the python interface, use `retrieve`
    instead"""
    tool = 'batch/'

    query = list(set(query.split('\n')))
    queries = [query[i:i+100] for i in range(0, len(query), 100)]

    data = {'format':format}

    responses = [requests.post(url + tool, data=data, files={'file':' '.join(query)}) for query in queries]
    page = ''.join([response.text for response in responses])
    return page

def retrieve(ids, format='txt'):
    """ request entries by uniprot acc using batch retrieval

    Args:
        query: list of ids to retrieve
        format: txt by default

    Help:
        possible formats:
        txt, xml, rdf, fasta, gff"""
    if type(ids) is not list:
        ids = [ids]
    return _retrieve(' '.join(ids), format)

def _map(query, f, t, format='tab'):
    """ _map is not meant for use with the python interface, use `map` instead
    """
    tool = 'mapping/'

    data = {
            'from':f,
            'to':t,
            'format':format,
            'query': query
            }
    response = requests.post(url + tool, data=data)
    page = response.text
    return page

def map(ids, f, t, format='tab'):
    """ map a list of ids from one format onto another using uniprots mapping api
    
    Args:
        query: id or list of ids to be mapped
        f: from ACC | P_ENTREZGENEID | ...
        t: to ...
        format: tab by default

    Help:
        for a list of all possible mappings visit
        'http://www.uniprot.org/faq/28'
    """
    if type(ids) is not list:
        ids = [ids]
    page = _map(' '.join(ids), f, t, format)
    result = dict()
    for row in page.splitlines()[1:]:
        key, value = row.split('\t')
        if key in result:
            result[key].add(value)
        else:
            result[key] = set([value])
    return result

if __name__ == '__main__':
    import argparse
    import sys

    parser = argparse.ArgumentParser(description='retrieve uniprot mapping')
    subparsers = parser.add_subparsers(dest='tool')

    mapping = subparsers.add_parser('map')
    mapping.add_argument('f', help='from')
    mapping.add_argument('t', help='to')
    mapping.add_argument('inp', nargs='?', type=argparse.FileType('r'),
            default=sys.stdin, help='input file (default: stdin)')
    mapping.add_argument('out', nargs='?', type=argparse.FileType('w'),
            default=sys.stdout, help='output file (default: stdout)')
    mapping.add_argument('--format', default='tab', help='output format')

    retrieve = subparsers.add_parser('retrieve')
    retrieve.add_argument('inp', metavar = 'in', nargs='?', type=argparse.FileType('r'),
            default=sys.stdin, help='input file (default: stdin)')
    retrieve.add_argument('out', nargs='?', type=argparse.FileType('w'),
            default=sys.stdout, help='output file (default: stdout)')
    retrieve.add_argument('-f', '--format', help='specify output format', default='txt')

    args = parser.parse_args()
    query = args.inp.read()

    if args.tool == 'map':
        args.out.write(_map(query, args.f, args.t, args.format))

    elif args.tool == 'retrieve':
        args.out.write(_retrieve(query, format=args.format))


