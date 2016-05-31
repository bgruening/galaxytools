"""Utils for PFAM."""

import xml.etree.ElementTree as ET
from math import ceil
from time import sleep
from xml.etree.ElementTree import ParseError

import requests

import pandas as pd

import fasta_utils

__author__ = "Gianluca Corrado"
__copyright__ = "Copyright 2016, Gianluca Corrado"
__license__ = "MIT"
__maintainer__ = "Gianluca Corrado"
__email__ = "gianluca.corrado@unitn.it"
__status__ = "Production"


def search_header():
    """Return the header of a Pfam scan search."""
    return "<seq id>        <alignment start>       <alignment end> \
    <envelope start>        <envelope end>  <hmm acc>       <hmm name>\
          <type>  <hmm start>     <hmm end>       <hmm length>    <bit score>\
               <E-value>       <significance>  <clan>\n"


def sequence_search(seq_id, seq):
    """
    Search a sequence against PFAM.

    Input
    -----
    seq_id : str
        Name of the protein sequence.
    seq : str
        Protein sequence.

    Output
    ------
    ret : str
        Formatted string containing the results of the Pfam scan for the
        given sequence
    """
    def add_spaces(text, mul=8):
        """Add spaces to a string."""
        l = len(text)
        next_mul = int(ceil(l / mul) + 1) * mul
        offset = next_mul - l
        if offset == 0:
            offset = 8
        return text + " " * offset

    url = "http://pfam.xfam.org/search/sequence"
    params = {'seq': seq,
              'evalue': '1.0',
              'output': 'xml'}
    req = requests.get(url, params=params)
    xml = req.text
    try:
        root = ET.fromstring(xml)
    # sometimes Pfam returns the HTML code
    except ParseError:
        print "resending: %s" % seq_id
        return "%s" % sequence_search(seq_id, seq)

    result_url = root[0][1].text
    # wait for Pfam to compute the results
    sleep(4)
    while True:
        req2 = requests.get(result_url)
        if req2.status_code == 200:
            break
        else:
            sleep(1)
    result_xml = req2.text
    root = ET.fromstring(result_xml)
    try:
        matches = root[0][0][0][0][:]
    # Sometimes raised when the sequence has no matches
    except IndexError:
        return ""
    ret = ""
    for match in matches:
        for location in match:
            ret += add_spaces(seq_id)
            ret += add_spaces(location.attrib['ali_start'])
            ret += add_spaces(location.attrib['ali_end'])
            ret += add_spaces(location.attrib['start'])
            ret += add_spaces(location.attrib['end'])
            ret += add_spaces(match.attrib['accession'])
            ret += add_spaces(match.attrib['id'])
            ret += add_spaces(match.attrib['class'])
            ret += add_spaces(location.attrib['hmm_start'])
            ret += add_spaces(location.attrib['hmm_end'])
            ret += add_spaces("None")
            ret += add_spaces(location.attrib['bitscore'])
            ret += add_spaces(location.attrib['evalue'])
            ret += add_spaces(location.attrib['significant'])
            ret += "None\n"
    return ret


def read_pfam_output(pfam_out_file):
    """Read the output of PFAM scan."""
    cols = ["seq_id", "alignment_start", "alignment_end", "envelope_start",
            "envelope_end", "hmm_acc", "hmm_name", "type", "hmm_start",
            "hmm_end", "hmm_length", "bit_score", "E-value", "significance",
            "clan"]
    try:
        data = pd.read_table(pfam_out_file,
                             sep="\s*", skip_blank_lines=True, skiprows=1,
                             names=cols, engine='python')
    except:
        return None
    return data


def download_seed_seqs(acc):
    """
    Download seed sequences from PFAM.

    Input
    -----
    acc : str
        Accession number of a Pfam domain

    Output
    ------
    fasta : str
        Seed sequences in fasta format
    """
    url = "http://pfam.xfam.org/family/%s/alignment/seed" % acc
    req = requests.get(url)
    stockholm = req.text
    fasta = fasta_utils.stockholm2fasta(stockholm)
    return fasta
