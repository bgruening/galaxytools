#!/usr/bin/env python

import io
import urllib2
import urllib
import httplib


def getListFromFile(infile):
    idlist = []
    for line in infile:
        line = line.strip()
        if line.isdigit():
            idlist.append(line)
    return idlist


def getListString(args):
    if args.id_type_ff == "file":
        list_string = ",".join(getListFromFile(open(args.id_value, "r")))
    else:
        list_string = args.id_value.strip().replace("__cr____cn__", ",")
    return list_string


def getresult(url):
    try:
        connection = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return ""
    else:
        return connection.read().rstrip()


def store_result_get(url, outfile):
    data = getresult(url)
    outfile.write(data)
    outfile.close()


def store_result_post(url, post, outfile):
    data = urllib.urlencode(post)
    headers = {"Content-Type": "application/x-www-form-urlencoded"}
    req = urllib2.Request(url, data, headers)
    response = urllib2.urlopen(req)
    the_page = response.read()
    outfile.write(the_page)
    outfile.close()
