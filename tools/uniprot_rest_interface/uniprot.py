#!/usr/bin/env python
"""
uniprot python interface
to access the uniprot database

Based on work from Jan Rudolph: https://github.com/jdrudolph/uniprot
available services:
    map
    retrieve
"""
import argparse
import sys

import requests

url = 'http://www.uniprot.org/'


def _retrieve(query, format='txt'):
    """_retrieve is not meant for use with the python interface, use `retrieve`
    instead"""
    tool = 'batch/'

    query = list(set(query.split('\n')))
    queries = [query[i:i+100] for i in range(0, len(query), 100)]

    data = {'format': format}

    responses = [requests.post(url + tool, data=data, files={'file': ' '.join(_)}) for _ in queries]
    page = ''.join(response.text for response in responses)
    return page


def _map(query, f, t, format='tab'):
    """ _map is not meant for use with the python interface, use `map` instead
    """
    tool = 'mapping/'

    data = {
            'from': f,
            'to': t,
            'format': format,
            'query': query
            }
    response = requests.post(url + tool, data=data)
    page = response.text
    return page


if __name__ == '__main__':
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
    retrieve.add_argument('inp', metavar='in', nargs='?', type=argparse.FileType('r'),
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
