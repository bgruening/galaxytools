#!/usr/bin/env python
"""
    Input: Molecules in SDF, SMILES ...
    Output: Moleculs filtered with specified substructures.
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

def parse_command_line():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--infile', required=True, help='Molecule file.')
    parser.add_argument('--iformat', help='Input format.')
    parser.add_argument('--fastsearch-index', dest="fastsearch_index", 
        required=True, help='Path to the openbabel fastsearch index.')
    parser.add_argument('-o', '--outfile', required=True, help='Path to the output file.')
    parser.add_argument('--oformat', 
        default='smi', help='Output file format')
    parser.add_argument("--max-candidates", dest="max_candidates", type=int,
                    default=4000, help="The maximum number of candidates.")
    parser.add_argument('-p', '--processors', type=int, 
        default=multiprocessing.cpu_count())
    return parser.parse_args()

results = list()
def mp_callback(res):
    results.append(res)

def mp_helper( query, args ):
    """
        Helper function for multiprocessing.
        That function is a wrapper around the following command:
        obabel file.fs -s"smarts" -Ooutfile.smi -al 999999999
    """

    if args.oformat == 'names':
        opts = '-osmi -xt'
    else:
        opts = '-o%s' % args.oformat

    tmp = tempfile.NamedTemporaryFile(delete=False)
    cmd = 'obabel -ifs %s -O %s %s -s%s -al %s' % (args.fastsearch_index, tmp.name, opts, query, args.max_candidates)

    child = subprocess.Popen(cmd.split(),
        stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = child.communicate()
    return_code = child.returncode

    if return_code:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        sys.stderr.write("Return error code %i from command:\n" % return_code)
        sys.stderr.write("%s\n" % cmd)
    else:
        sys.stdout.write(stdout)
        sys.stdout.write(stderr)
    return (tmp.name, query)


def get_smiles_or_smarts( args ):
    """
    Wrapper to retrieve a striped SMILES or SMARTS string from different input formats.
    """
    if args.iformat in ['smi', 'text', 'tabular']:
        with open( args.infile ) as text_file:
            for line in text_file:
                yield line.split('\t')[0].strip()
    else:
        # inchi or sdf files
        for mol in pybel.readfile( args.iformat, args.infile ):
            yield mol.write('smiles').split('\t')[0]

def substructure_search( args ):

    pool = multiprocessing.Pool( args.processors )
    for query in get_smiles_or_smarts( args ):
        pool.apply_async(mp_helper, args=(query, args), callback=mp_callback)
        #mp_callback( mp_helper(query, args) )
    pool.close()
    pool.join()

    if args.oformat == 'names':
        out_handle = open( args.outfile, 'w' )
        for result_file, query in results:
            with open(result_file) as res_handle:
                for line in res_handle:
                    out_handle.write('%s\t%s\n' % ( line.strip(), query ))
            os.remove( result_file )
        out_handle.close()
    else:
        out_handle = open( args.outfile, 'wb' )
        for result_file, query in results:
            res_handle = open(result_file,'rb')
            shutil.copyfileobj( res_handle, out_handle )
            res_handle.close()
            os.remove( result_file )
        out_handle.close()


def __main__():
    """
        Multiprocessing Open Babel Substructure Search.
    """
    args = parse_command_line()
    substructure_search( args )

if __name__ == "__main__" :
    __main__()
