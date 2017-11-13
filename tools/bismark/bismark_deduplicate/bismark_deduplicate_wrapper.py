#!/usr/bin/python

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile
import logging
from glob import glob

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tool_dir', dest='tool_dir', action='store', nargs=1, metavar='tool_dir', type=str)
    parser.add_argument('-p', action='store_true')
    parser.add_argument('-s', action='store_true')
    parser.add_argument('--input', dest='input', action='store', nargs=1, metavar='input', type=str)
    parser.add_argument('--output_report', dest='output_report', action='store', nargs=1, metavar='output_report', type=str)
    parser.add_argument('--output_bam', dest='output_bam', action='store', nargs=1, metavar='output_report', type=str)
    parser.add_argument('--log_report', dest='log_report', action='store', nargs=1, metavar='log_filename', type=str)
    args = parser.parse_args()
    return args

def __main__():
    args = get_arg()

    tmp_dir = tempfile.mkdtemp(prefix='tmp', suffix='')
    os.chdir(tmp_dir)

    if args.log_report:
        logging.basicConfig(level=logging.INFO, filename=args.log_report[0], filemode="a+", format='%(message)s')
    else:
        logging.basicConfig(level=logging.INFO, filename=os.path.join(tmp_dir, 'log_report.txt'), filemode="a+", format='%(message)s')

    default_reads_name = 'submitted_reads.bam'
    os.symlink(args.input[0], default_reads_name)

    if args.p is True:
        sPaired = '-p'
    if args.s is True:
        sPaired = '-s'

    cmd = 'perl %s %s duplicated_reads.bam --bam' % (os.path.join(args.tool_dir[0], 'deduplicate_bismark'), sPaired)
    logging.info('COMMAND LINE:\n\n%s' % cmd)

    proc = subprocess.Popen(['perl', os.path.join(args.tool_dir[0], 'deduplicate_bismark'), sPaired, default_reads_name, '--bam'], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    proc_out, proc_err = proc.communicate()

    logging.info("__________________________________________________________________\n")
    logging.info("BISMARK DEDUPLICATE STDOUT:\n\n%s" % proc_out)
    if proc_err:
        logging.critical("__________________________________________________________________\n")
        logging.critical("BISMARK DEDUPLICATE WARNING:\n\n%s" % proc_err)
        sys.exit("Dedpulicate Bismark crashed with the folowing error message:\n%s" % proc_err)

    shutil.move( glob('*deduplicated.bam')[0], args.output_bam[0] )
    shutil.move( glob('*deduplication_report.txt')[0], args.output_report[0])

    cleanup_before_exit(tmp_dir)

if __name__=="__main__": __main__()