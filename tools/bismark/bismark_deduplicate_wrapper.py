#!/usr/bin/python

import argparse
import logging
import os
import shutil
import subprocess
import sys
import signal
import tempfile
from glob import glob


def stop_err(logger, msg):
    logger.critical(msg)
    sys.exit(1)


def restore_sigpipe():
    """
    Needed to handle samtools view
    """
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)


def log_subprocess_output(logger, pipe):
    for line in iter(pipe.readline, b''):
        logger.debug(line.decode().rstrip())


def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--single_or_paired',  dest='single_or_paired')
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

    # ensure the input has a .bam suffix
    tmp_dir = tempfile.mkdtemp(prefix='tmp', suffix='')
    os.chdir(tmp_dir)
    default_reads_name = 'submitted_reads.bam'
    os.symlink(args.input, default_reads_name)

    single_or_paired = '-s' if args.single_or_paired == 'single' else '-p'
    cmd = ['deduplicate_bismark', single_or_paired, default_reads_name, '--bam']
    logger.info("Deduplicating with: '%s'", " ".join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               preexec_fn=restore_sigpipe)
    proc_out, proc_err = process.communicate()
    logger.info(proc_out)
    if process.returncode != 0:
        stop_err(logger, "Bismark deduplication error (also check the log file if any)!\n%s" % proc_err)

    deduplicated_out_name = 'submitted_reads.deduplicated.bam'
    deduplicated_report_name = 'submitted_reads.deduplication_report.txt'
    logger.debug("Moving '%s' to galaxy: '%s'.", deduplicated_out_name, args.output_bam)
    shutil.move(deduplicated_out_name, args.output_bam )
    logger.debug("Moving '%s' to galaxy: '%s'.", deduplicated_report_name, args.output_report)
    shutil.move('submitted_reads.deduplication_report.txt', args.output_report)
    logger.debug("Done.")

if __name__=="__main__": __main__()