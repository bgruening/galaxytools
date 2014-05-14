#!/usr/bin/env python


import sys, os
import tempfile
import argparse
import urllib2, urllib, httplib
import readfile
import xml.sax as sax

#parse a file where there is some AID/CID/SID as key and several other IDs belonging
#to it
class DictHandler(sax.handler.ContentHandler):

    def __init__(self): 
        self.ergebnis = {} 
        self.schluessel = "" 
        self.wert = "" 
        self.aktiv = None
        self.keywert = "AID"
        self.valuewert ="CID"

    def __init__(self, key, value): 
        self.ergebnis = {} 
        self.schluessel = "" 
        self.wert = "" 
        self.aktiv = None
        self.keywert = key
        self.valuewert = value

    def startElement(self, name, attrs): 
        if name == "Information": 
            self.schluessel = "" 
            self.wert = "" 
        elif name == self.keywert or name== self.valuewert: 
            self.aktiv = name 

    def endElement(self, name): 
        if name == self.keywert:
            self.schluessel=self.schluessel.strip()
            self.ergebnis[self.schluessel]=[]
            #print("huhn")
            self.aktiv=None
        elif name == self.valuewert: 
            self.aktiv = None
            self.ergebnis[self.schluessel].append(self.wert)
            self.wert=""
    def characters(self, content): 
        if self.aktiv == self.keywert: 
            self.schluessel += content 
        elif self.aktiv == self.valuewert:
            self.wert += content
            

def dict_from_xml(xmlfile, key, value):
    handler = DictHandler(key, value) 
    parser = sax.make_parser() 
    parser.setContentHandler(handler)
    parser.parse(xmlfile)
    dic=handler.ergebnis
    return dic
    
#get every self.keywert as a list
#returns a dictionary with self.keywert as key and as value the list of self.valuewerts
def getAllAssayIDs():
    url="http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/type/all/aids/TXT"
    data=readfile.getresult(url)
    self.keywertlist=readfile.getListFromString(data)
    return self.keywertlist

#get all ids in the url xml file and produce a dictionary returned by the xml parser
def get_dict_key_value(url, key, value):
    xml=readfile.getresult(url)
    tmp = tempfile.TemporaryFile() 
    tmp.write(xml)
    tmp.seek(0)
    dic=dict_from_xml(tmp, key, value)
    tmp.close()
    return dic

def getIDofLine(line):
    arr=line.split(">")
    if len(arr) > 1:
        self.keywert=arr[1].split("<")[0]
        return self.keywert
    else:
        return "-1"

def write_to_sf(iddict, outfile, sep):
    for key_id in iddict:
        for val_id in iddict[key_id]:
            outfile.write(key_id)
            outfile.write(sep)
            outfile.write(val_id)
            outfile.write("\n")

