#!/usr/bin/env python


import sys, os
import argparse
import urllib2, urllib, httplib
import readfile
import xml.sax as sax

class DictHandler(sax.handler.ContentHandler):

    def __init__(self): 
        self.ergebnis = {} 
        self.schluessel = "" 
        self.wert = "" 
        self.aktiv = None 

    def startElement(self, name, attrs): 
        if name == "Information": 
            self.schluessel = "" 
            self.wert = "" 
        elif name == "AID" or name=="CID": 
            self.aktiv = name 

    def endElement(self, name): 
        if name == "AID":
            self.schluessel=self.schluessel.strip()
            self.ergebnis[self.schluessel]=[]
            #print("huhn")
            self.aktiv=None
        elif name == "CID": 
            self.aktiv = None
            self.ergebnis[self.schluessel].append(self.wert)
            self.wert=""
    def characters(self, content): 
        if self.aktiv == "AID": 
            self.schluessel += content 
        elif self.aktiv == "CID":
            self.wert += content
            

def give_aid_cid_dict_from_xml(xmlfile):
    handler = DictHandler() 
    parser = sax.make_parser() 
    parser.setContentHandler(handler)
    parser.parse(xmlfile)
    dic=handler.ergebnis
    return dic
    
#get every aid as a list
#returns a dictionary with aid as key and as value the list of cids
def getAllAssayIDs():
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/type/all/aids/TXT"
    data=readfile.getresult(url)
    aidlist=readfile.getListFromString(data)
    return aidlist


def getIDofLine(line):
    arr=line.split(">")
    if len(arr) > 1:
        aid=arr[1].split("<")[0]
        return aid
    else:
        return "-1"

def write_to_csv(aid_cid_dict, outfile):
    for key in aid_cid_dict:
        for cid in aid_cid_dict[key]:
            outfile.write(key)
            outfile.write(",")
            outfile.write(cid)
            outfile.write("\n")

