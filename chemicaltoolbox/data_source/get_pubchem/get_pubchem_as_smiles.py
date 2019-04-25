#!/usr/bin/env python

__author__ = 'Bjoern Gruening'
__version__ = '0.1'
__date__ = '2012'
__license__ = 'GLP3+'

import ftplib
import os, sys
import argparse
import subprocess
from multiprocessing import Pool
import tempfile
import shutil


def main(output, processors = 4):
    output_handle = open(output,'w+')

    td = tempfile.mkdtemp()
    ftp = ftplib.FTP('ftp.ncbi.nih.gov')
    ftp.login()
    ftp.cwd('/pubchem/Compound/CURRENT-Full/SDF/')
    filelist = ftp.nlst()
    
    pool = Pool(processes = processors)
    filenames = zip(filelist, [td]*len(filelist))
    result = pool.map_async(fetch_convert, filenames)
    result.get()

    for filename in os.listdir(td):
        path = os.path.join(td, filename)
        shutil.copyfileobj(open(path, 'rb'), output_handle)

    output_handle.close()
    shutil.rmtree(td)

def fetch_convert(args):
    (filename, td) = args
    tmp_name = os.path.join( td, filename)
    subprocess.call( ['wget', '-O', tmp_name, os.path.join('ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/SDF/', filename)] )
    output = os.path.join(td, filename) + '.smi'
    subprocess.call(["obabel", "-isdf", tmp_name, "-ocan", '-O', output])
    os.remove(tmp_name)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Download the whole PubChem and converts it to canonical SMILES on the fly.')
    parser.add_argument("-o", "--output", dest="output",
                    required=True,
                    help="Path to the output file.")
    parser.add_argument("-p", "--processors", dest="processors",
                    type=int, default=10,
                    help="How many processors you want to use.")

    options = parser.parse_args()
    main( options.output, options.processors )

