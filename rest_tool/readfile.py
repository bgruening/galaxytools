#!/usr/bin/env python

import io
import urllib2, urllib, httplib
def getListFromFile(file):
    idlist=[]
    for line in file:
        if int(line):
            idlist.append(line.strip())
    return idlist

def getresult(url):
    try:
        connection = urllib2.urlopen(url)
    except urllib2.HTTPError, e:
        return ""
    else:
        return connection.read().rstrip()
        
