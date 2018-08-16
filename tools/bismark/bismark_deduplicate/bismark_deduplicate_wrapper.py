#!/usr/bin/python

import argparse
import logging
import os
import shutil
import subprocess
import sys
import tempfile
from glob import glob


def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)


def stop_err(logger, msg):
    logger.critical(msg)
    sys.exit(1)


def log_subprocess_output(logger, pipe):
    for line in iter(pipe.readline, b''):
        logger.debug(line.decode().rstrip())


def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action='store_true')
    parser.add_argument('-s', action='store_true')
    parser.add_argument('--input', dest='input', metavar='input')
    parser.add_argument('--output_report', dest='output_report', metavar='output_report')
    parser.add_argument('--output_bam', dest='output_bam', metavar='output_report')
    parser.add_argument('--log_report', dest='log_report', metavar='log_filename', type=str)
    args = parser.parse_args()
    return args


def __main__():
    args = get_arg()

    logger = logging.getLogger('bismark_deduplicate_wrapper')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    if args.log_report:
        ch.setLevel(logging.WARNING)
        handler = logging.FileHandler(args.log_report)
        handler.setLevel(logging.DEBUG)
        logger.addHandler(handler)
    else:
        ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    tmp_dir = tempfile.mkdtemp(prefix='tmp', suffix='')
    os.chdir(tmp_dir)

    default_reads_name = 'submitted_reads.bam'
    os.symlink(args.input, default_reads_name)

    if args.p is True:
        sPaired = '-p'
    if args.s is True:
        sPaired = '-s'

    cmd = ['deduplicate_bismark', sPaired, default_reads_name, '--bam']
    logger.info("Deduplicating with: '%s'", " ".join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(logger, process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        stop_err(logger, "Bismark deduplication error (also check the log file if any)!\n%s" % process.stderr)

    shutil.move( glob('*deduplicated.bam')[0], args.output_bam )
    shutil.move( glob('*deduplication_report.txt')[0], args.output_report)

    cleanup_before_exit(tmp_dir)

if __name__=="__main__": __main__()