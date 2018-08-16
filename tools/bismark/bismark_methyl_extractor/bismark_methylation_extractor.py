#!/usr/bin/env python

import argparse
import logging
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile
from glob import glob


def stop_err(logger, msg):
    logger.critical(msg)
    sys.exit(1)


def log_subprocess_output(logger, pipe):
    for line in iter(pipe.readline, b''):
        logger.debug(line.decode().rstrip())


def zipper(dir, zip_file):
    output_files_regex = re.compile('^(Non_)?C[pH][GH]_.*')
    bedgraph_regex = re.compile('.*bedGraph.gz')
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        archive_root = os.path.abspath(root)[root_len:]
        for f in files:
            if re.search(output_files_regex, f) or re.search(bedgraph_regex, f):
                fullpath = os.path.join(root, f)
                archive_name = os.path.join(archive_root, f)
                zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


def build_genome_dir(genome_file):
    tmp_genome_dir = tempfile.mkdtemp(prefix='tmp')
    genome_path = os.path.join(tmp_genome_dir, '.'.join(os.path.split(genome_file)[1].split('.')[:-1]))
    try:
        # Create a hard link pointing to genome_file named 'genome_path'.fa.
        os.symlink(genome_file, genome_path + '.fa')
    except Exception as e:
        if os.path.exists(tmp_genome_dir):
            shutil.rmtree(tmp_genome_dir)
        stop_err('Error in linking the reference database!\n%s' % e)
    return tmp_genome_dir


def __main__():
    # Parse Command Line
    parser = argparse.ArgumentParser(description='Wrapper for the bismark methylation caller.')

    # input options
    parser.add_argument('--infile', help='Input file in SAM or BAM format.')
    parser.add_argument('--single-end', dest='single_end', action="store_true")
    parser.add_argument('--paired-end', dest='paired_end', action="store_true")

    parser.add_argument('--multicore', dest='multicore', type=int, default=1)
    parser.add_argument('--splitting_report', dest='splitting_report')
    parser.add_argument('--mbias_report', dest='mbias_report')
    parser.add_argument('--cytosine_report', dest="cytosine_report")
    parser.add_argument('--genome_file', dest="genome_file")
    parser.add_argument('--cx_context', action="store_true")

    parser.add_argument('--comprehensive', action="store_true")
    parser.add_argument('--merge-non-cpg', dest='merge_non_cpg', action="store_true")
    parser.add_argument('--no-overlap', dest='no_overlap', action="store_true")
    parser.add_argument('--compress')
    parser.add_argument('--ignore', dest='ignore', type=int)
    parser.add_argument('--ignore_r2', dest='ignore_r2', type=int)
    parser.add_argument('--ignore_3prime', dest='ignore_3prime', type=int)
    parser.add_argument('--ignore_3prime_r2', dest='ignore_3prime_r2', type=int)
    parser.add_argument('--log_report', dest='log_report', metavar='log_filename', type=str)
    args = parser.parse_args()

    logger = logging.getLogger('bismark_methylation_extractor_wrapper')
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

    # Build methylation extractor command
    output_dir = tempfile.mkdtemp()
    cmd = ['bismark_methylation_extractor', '--multicore', str(args.multicore), '--no_header', '-o', output_dir]

    # Set up all options
    if args.single_end:
        cmd.append('--single-end')
    else:
        cmd.append('--paired-end')
    if args.no_overlap:
        cmd.append('--no_overlap')
    if args.ignore:
        cmd.extend(['--ignore', str(args.ignore)])
    if args.ignore_r2:
        cmd.extend(['--ignore_r2', str(args.ignore_r2)])
    if args.ignore_3prime:
        cmd.extend(['--ignore_3prime', str(args.ignore_3prime)])
    if args.ignore_3prime_r2:
        cmd.extend(['--ignore_3prime_r2', str(args.ignore_3prime_r2)])
    if args.comprehensive:
        cmd.append('--comprehensive')
    if args.merge_non_cpg:
        cmd.append('--merge_non_CpG')
    if args.splitting_report:
        cmd.append('--report')
    tmp_genome_dir = None
    if args.cytosine_report:
        tmp_genome_dir = build_genome_dir(args.genome_file)
        if args.cx_context:
            cmd.extend(
                ['--bedgraph', '--CX_context', '--cytosine_report', '--CX_context', '--genome_folder', tmp_genome_dir])
        else:
            cmd.extend(['--bedgraph', '--cytosine_report', '--genome_folder', tmp_genome_dir])

    cmd.append(args.infile)

    # Run
    logger.info("Methylation extractor run with: '%s'", " ".join(cmd))
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    with process.stdout:
        log_subprocess_output(logger, process.stdout)
    exitcode = process.wait()
    if exitcode != 0:
        stop_err(logger, "Bismark methylation extractor error (also check the log file if any)!\n%s" % process.stderr)

    # collect and copy output files
    zipper(output_dir, args.compress)

    # cytosine report
    if args.cytosine_report:
        if args.cx_context:
            shutil.move(glob(os.path.join(output_dir, '*CX_report.txt'))[0], args.cytosine_report)
        else:
            shutil.move(glob(os.path.join(output_dir, '*CpG_report.txt'))[0], args.cytosine_report)
    # splitting report
    if args.splitting_report:
        shutil.move(glob(os.path.join(output_dir, '*_splitting_report.txt'))[0], args.splitting_report)
    if args.mbias_report:
        shutil.move(glob(os.path.join(output_dir, '*M-bias.txt'))[0], args.mbias_report)

    # Clean up temp dirs
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    if tmp_genome_dir and os.path.exists(tmp_genome_dir):
        shutil.rmtree(tmp_genome_dir)


if __name__ == "__main__": __main__()
