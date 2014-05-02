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
        
def store_result_get(url, outfile):
    data=getresult(url)
    outfile.write(data)
    outfile.close()

def store_result_post(url, post, outfile):
    data = urllib.urlencode(post)
    headers={"Content-Type" : "application/x-www-form-urlencoded"}
    req = urllib2.Request(url, data, headers)
    response = urllib2.urlopen(req)
    the_page = response.read()
    outfile.write(the_page)
    outfile.close()
