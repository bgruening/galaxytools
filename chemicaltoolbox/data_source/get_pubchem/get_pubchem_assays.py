#!/usr/bin/env python

__author__ = 'Bjoern Gruening'
__version__ = '0.1'
__date__ = '2014'
__license__ = 'GLP3+'

import ftplib
import os, sys
import argparse
import subprocess
from multiprocessing import Pool
import tempfile
import shutil
import urllib
import zipfile
import gzip


def main(output, processors = 4, white_list = ['Active','Inconclusive', 'Inactive']):
    """
        Starting multiple processes to download and extract PubChem Assay data.
    """
    td = tempfile.mkdtemp()
    ftp = ftplib.FTP('ftp.ncbi.nih.gov')
    ftp.login()
    ftp.cwd('/pubchem/Bioassay/Concise/CSV/Data/')
    filelist = ftp.nlst()

    pool = Pool(processes = processors)
    triplestore = zip(filelist, [td]*len(filelist), [white_list]*len(filelist))

    result = pool.map_async(fetch_convert, triplestore)
    result.get()

    with open(output,'w+') as output_handle:
        for filename in os.listdir( td ):
            path = os.path.join( td, filename )
            shutil.copyfileobj(open(path, 'rb'), output_handle)

    shutil.rmtree( td )

def fetch_convert(args):
    (filename, td, white_list) = args
    tmp_name = os.path.join( td, filename)
    urllib.urlretrieve(os.path.join('ftp://ftp.ncbi.nih.gov/pubchem/Bioassay/Concise/CSV/Data/', filename), tmp_name)

    temp_dir = tempfile.mkdtemp()
    with zipfile.ZipFile(tmp_name, "r") as z:
        z.extractall(temp_dir)

    output = os.path.join(td, filename) + '.tsv'
    with open(output, 'w+') as out_handle:
        for root, dirs, files in os.walk( temp_dir ):
            for filename in files:
                gzfile_path = os.path.join( root, filename )
                with gzip.open(gzfile_path, 'rb') as gzfile:
                    gzfile.readline() # skip first line
                    for line in gzfile:
                        cols = line.split(',')
                        PUBCHEM_ACTIVITY_OUTCOME = cols[2]
                        if PUBCHEM_ACTIVITY_OUTCOME in white_list:
                            out_handle.write( '%s' % line.replace(',', '\t') )

    os.remove(tmp_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download the whole PubChem and converts it to canonical SMILES on the fly.')
    parser.add_argument("-o", "--output", dest="output",
                    required=True,
                    help="Path to the output file.")
    parser.add_argument("-p", "--processors", dest="processors",
                    type=int, default=10,
                    help="How many processors you want to use.")
    parser.add_argument("-w", "--white-list", dest="white_list",
                    default="Active,Inconclusive,Inactive",
                    help="List of comma separated PUBCHEM_ACTIVITY_OUTCOME values that should be fetched.")

    options = parser.parse_args()
    main( options.output, options.processors, options.white_list.split(',') )

