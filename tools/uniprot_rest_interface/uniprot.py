#!/usr/bin/env python
"""
uniprot python interface
to access the uniprot database

Based on work from Jan Rudolph: https://github.com/jdrudolph/uniprot
available services:
    map
    retrieve

rewitten using inspiration form: https://findwork.dev/blog/advanced-usage-python-requests-timeouts-retries-hooks/
"""
import argparse
import sys

import requests
from requests.adapters import HTTPAdapter
from requests.packages.urllib3.util.retry import Retry


DEFAULT_TIMEOUT = 5  # seconds
URL = 'https://www.uniprot.org/'

retry_strategy = Retry(
    total=5,
    backoff_factor=2,
    status_forcelist=[429, 500, 502, 503, 504],
    allowed_methods=["HEAD", "GET", "OPTIONS", "POST"]
)


class TimeoutHTTPAdapter(HTTPAdapter):
    def __init__(self, *args, **kwargs):
        self.timeout = DEFAULT_TIMEOUT
        if "timeout" in kwargs:
            self.timeout = kwargs["timeout"]
            del kwargs["timeout"]
        super().__init__(*args, **kwargs)

    def send(self, request, **kwargs):
        timeout = kwargs.get("timeout")
        if timeout is None:
            kwargs["timeout"] = self.timeout
        return super().send(request, **kwargs)


def _map(query, f, t, format='tab', chunk_size=100):
    """ _map is not meant for use with the python interface, use `map` instead
    """
    tool = 'uploadlists/'
    data = {'format': format, 'from': f, 'to': t}

    req = []
    for i in range(0, len(query), chunk_size):
        q = query[i:i + chunk_size]
        req.append(dict([("url", URL + tool),
                         ('data', data),
                         ("files", {'file': ' '.join(q)})]))
    return req
    response = requests.post(URL + tool, data=data)
    response.raise_for_status()
    page = response.text
    if "The service is temporarily unavailable" in page:
        exit("The UNIPROT service is temporarily unavailable. Please try again later.")
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

    # get the IDs from the file as sorted list
    # (sorted is convenient for testing)
    query = set()
    for line in args.inp:
        query.add(line.strip())
    query = sorted(query)

    if args.tool == 'map':
        pload = _map(query, args.f, args.t, chunk_size=100)
    elif args.tool == 'retrieve':
        pload = _map(query, 'ACC+ID', 'ACC', args.format, chunk_size=100)

    adapter = TimeoutHTTPAdapter(max_retries=retry_strategy)
    http = requests.Session()
    http.mount("https://", adapter)
    for i, p in enumerate(pload):
        response = http.post(**p)
        args.out.write(response.text)
    http.close()