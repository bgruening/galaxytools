#!/usr/bin/env python
"""
Description:
    Executes the tbl2asn tool.

Usage:
    tbl2asn.py -f <file_type> -o <output_type> -t <tbl_file> -a <fsa_file> -s <sbt_file>

Reference:
    http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/

"""
import sys
import os
import argparse
import re


def main(filetype, output, tblFile, fsaFile, sbtFile):
    os.system("mkdir ncbi_tbl_temp")
    os.system("cp " + tblFile + " ./ncbi_tbl_temp/temp_file_tbl2asn.tbl")
    os.system("cp " + fsaFile + " ./ncbi_tbl_temp/temp_file_tbl2asn.fsa")
    outOpt = ""
    for i in output:
        outOpt +=  i
    if not sbtFile:
        os.system("tbl2asn -p ./ncbi_tbl_temp -a " + filetype + " -V " + outOpt)
    else:
        os.system("tbl2asn -t " + sbtFile + " -p ./ncbi_tbl_temp -a " + filetype + " -V " + outOpt)


if __name__ == '__main__':

    parser = argparse.ArgumentParser(prog = 'tbl2asn.py', description='Executes the tbl2asn tool.', prefix_chars='-+', epilog="")
    parser.add_argument('--version', action='version', version='%(prog)s 0.1')
    parser.add_argument('--filetype', '-f', dest='filetype', required=True, help='specifies file type')
    parser.add_argument('--output', '-o', dest='output', required=True, help='specifies output files (validation)')
    parser.add_argument('--tblFile', '-t', dest='tblFile', required=True, help='specifies the tbl file')
    parser.add_argument('--fsaFile', '-a', dest='fsaFile', required=True, help='specifies the fsa file')
    parser.add_argument('--sbtFile', '-s', dest='sbtFile', required=False, help='specifies sbt file')

    options = parser.parse_args()
    main(options.filetype, options.output, options.tblFile, options.fsaFile, options.sbtFile)
