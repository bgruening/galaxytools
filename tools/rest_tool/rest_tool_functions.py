#!/usr/bin/env python

import tempfile
import xml.sax as sax

import httplib
import readfile
import urllib2


class DictHandler(sax.handler.ContentHandler):
    """
    Parse XML files with AID/CID/SID as key and mappings to several other IDs
    """

    def __init__(self):
        self.ergebnis = {}
        self.schluessel = ""
        self.wert = ""
        self.aktiv = None
        self.keywert = "AID"
        self.valuewert = "CID"

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
        elif name == self.keywert or name == self.valuewert:
            self.aktiv = name

    def endElement(self, name):
        if name == self.keywert:
            self.schluessel = self.schluessel.strip()
            self.ergebnis[self.schluessel] = []
            self.aktiv = None
        elif name == self.valuewert:
            self.aktiv = None
            self.ergebnis[self.schluessel].append(self.wert)
            self.wert = ""

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
    dic = handler.ergebnis
    return dic


def getAllAssayIDs():
    url = "http://pubchem.ncbi.nlm.nih.gov/rest/pug/assay/type/all/aids/TXT"
    data = readfile.getresult(url)
    aids = readfile.getListFromString(data)
    return aids


def get_dict_key_value(url, key, value):
    """
    Get all ids from a XML file, obtained from the given URL
    and convert the information into a dictionary
    """
    xml = readfile.getresult(url)
    tmp = tempfile.TemporaryFile()
    tmp.write(xml)
    tmp.seek(0)
    dic = dict_from_xml(tmp, key, value)
    tmp.close()
    return dic


def getIDofLine(line):
    arr = line.split(">")
    if len(arr) > 1:
        self.keywert = arr[1].split("<")[0]
        return self.keywert
    else:
        return "-1"


def write_to_sf(iddict, outfile, sep):
    for key_id, values in iddict:
        for val_id in values:
            outfile.write("%s%s%s\n" % (key_id, sep, val_id))
