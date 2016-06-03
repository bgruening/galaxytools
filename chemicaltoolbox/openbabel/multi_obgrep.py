#!/usr/bin/env python
"""
    Input: Molecules in SDF, SMILES ...
    Output: Molecule file filtered with obgrep.
    Copyright 2013, Bjoern Gruening and Xavier Lucas
"""
import sys, os
import argparse
import openbabel
openbabel.obErrorLog.StopLogging()
import pybel
import multiprocessing
import tempfile
import subprocess
import shutil
import shlex

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Molecule file.')
    parser.add_argument('-q', '--query',  required=True, help='Query file, containing different SMARTS in each line.')
    parser.add_argument('-o', '--outfile', required=True, help='Path to the output file.')
    parser.add_argument("--iformat", help="Input format, like smi, sdf, inchi")
    parser.add_argument("--n-times", dest="n_times", type=int,
                    default=0, help="Print a molecule only if the pattern occurs # times inside the molecule.")
    parser.add_argument('-p', '--processors', type=int, default=multiprocessing.cpu_count())
    parser.add_argument("--invert-matches", dest="invert_matches", action="store_true",
                    default=False, help="Invert the matching, print non-matching molecules.")
    parser.add_argument("--only-name", dest="only_name", action="store_true",
                    default=False, help="Only print the name of the molecules.")
    parser.add_argument("--full-match", dest="full_match", action="store_true",
                    default=False, help="Full match, print matching-molecules only when the number of heavy atoms is also equal to the number of atoms in the SMARTS pattern.")
    parser.add_argument("--number-of-matches", dest="number_of_matches", action="store_true",
                    default=False, help="Print the number of matches.")
    return parser.parse_args()

results = list()
def mp_callback(res):
    results.append(res)

def mp_helper( query, args ):
    """
        Helper function for multiprocessing.
        That function is a wrapper around obgrep.
    """

    cmd_list = []
    if args.invert_matches:
        cmd_list.append('-v')
    if args.only_name:
        cmd_list.append('-n')
    if args.full_match:
        cmd_list.append('-f')
    if args.number_of_matches:
        cmd_list.append('-c')
    if args.n_times:
        cmd_list.append('-t %s' % str(args.n_times))

    tmp = tempfile.NamedTemporaryFile(delete=False)
    cmd = 'obgrep %s "%s" %s' % (' '.join(cmd_list), query, args.infile)
    child = subprocess.Popen(shlex.split(cmd),
        stdout=open(tmp.name, 'w+'), stderr=subprocess.PIPE)

    stdout, stderr = child.communicate()
    return (tmp.name, query)


def obgrep( args ):

    temp_file = tempfile.NamedTemporaryFile()
    temp_link = "%s.%s" % (temp_file.name, args.iformat)
    temp_file.close()
    os.symlink(args.infile, temp_link)
    args.infile = temp_link

    pool = multiprocessing.Pool( args.processors )
    for query in open( args.query ):
        pool.apply_async(mp_helper, args=(query.strip(), args), callback=mp_callback)
        #mp_callback( mp_helper(query.strip(), args) )
    pool.close()
    pool.join()

    out_handle = open( args.outfile, 'wb' )
    for result_file, query in results:
        res_handle = open(result_file,'rb')
        shutil.copyfileobj( res_handle, out_handle )
        res_handle.close()
        os.remove( result_file )
    out_handle.close()

    os.remove( temp_link )

def __main__():
    """
        Multiprocessing obgrep search.
    """
    args = parse_command_line()
    obgrep( args )

if __name__ == "__main__" :
    __main__()
