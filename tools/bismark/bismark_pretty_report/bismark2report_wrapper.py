#!/usr/bin/python

import argparse
import logging
import subprocess
import sys


def stop_err(logger, msg):
    logger.critical(msg)
    sys.exit()


def log_subprocess_output(logger, pipe):
    for line in iter(pipe.readline, b''):
        logger.debug(line.decode().rstrip())


def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--alignment_report', dest='alignment_report', action='store', metavar='alignment_report',
                        type=str)
    parser.add_argument('--dedup_report', dest='dedup_report', action='store', metavar='dedup_report', type=str)
    parser.add_argument('--splitting_report', dest='splitting_report', action='store', metavar='splitting_report',
                        type=str)
    parser.add_argument('--mbias_report', dest='mbias_report', action='store', metavar='mbias_report', type=str)
    parser.add_argument('--nucleotide_report', dest='nucleotide_report', action='store', metavar='nucleotide_report',
                        type=str)
    parser.add_argument('--output_html_report', dest='output_html_report', action='store', metavar='output_html_report',
                        type=str)
    parser.add_argument('--log_report', dest='log_report', action='store', metavar='log_report', type=str)
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

    alignment_option = '--alignment_report'
    alignment_report = args.alignment_report
    if args.dedup_report:
        dedup_option = '--dedup_report'
        dedup_report = args.dedup_report
    else:
        dedup_option = ''
        dedup_report = ''
    if args.splitting_report:
        splitting_option = '--splitting_report'
        splitting_report = args.splitting_report
    else:
        splitting_option = ''
        splitting_report = ''
    if args.mbias_report:
        mbias_option = '--mbias_report'
        mbias_report = args.mbias_report
    else:
        mbias_option = ''
        mbias_report = ''
    if args.nucleotide_report:
        nucleotide_option = '--nucleotide_report'
        nucleotide_report = args.nucleotide_report
    else:
        nucleotide_option = ''
        nucleotide_report = ''

    cmd = ['bismark2report', alignment_option, alignment_report, dedup_option, dedup_report, splitting_option,
           splitting_report, mbias_option, mbias_report, nucleotide_option, nucleotide_report,
           '--output', args.output_html_report]
    logger.info("Generating report with: '%s'", " ".join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(logger, process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        stop_err("Bismark pretty report error (also check the log file if any)!\n%s" % process.stderr)


if __name__ == "__main__": __main__()
