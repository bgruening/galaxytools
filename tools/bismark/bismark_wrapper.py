#!/usr/bin/env python

import argparse
import fileinput
import logging
import math
import os
import shutil
import subprocess
import sys
import tempfile
from glob import glob


def stop_err(logger, msg):
    logger.critical(msg)
    sys.exit(1)


def log_subprocess_output(logger, pipe):
    for line in iter(pipe.readline, b''):
        logger.debug(line.decode().rstrip())


def __main__():
    # Parse Command Line
    parser = argparse.ArgumentParser(description='Wrapper for the bismark bisulfite mapper.')
    parser.add_argument('-p', '--num-threads', dest='num_threads',
                        type=int, default=4, help='Use this many threads to align reads. The default is 4.')

    # input options
    parser.add_argument('--own-file', dest='own_file', help='')
    parser.add_argument('-D', '--indexes-path', dest='index_path',
                        help='Indexes directory; location of .ebwt and .fa files.')
    parser.add_argument('-O', '--output', dest='output')

    parser.add_argument('--output-report-file', dest='output_report_file')
    parser.add_argument('--suppress-header', dest='suppress_header', action="store_true")

    parser.add_argument('--mate-paired', dest='mate_paired', action='store_true', help='Reads are mate-paired',
                        default=False)

    parser.add_argument('-1', '--mate1', dest='mate1',
                        help='The forward reads file in Sanger FASTQ or FASTA format.')
    parser.add_argument('-2', '--mate2', dest='mate2',
                        help='The reverse reads file in Sanger FASTQ or FASTA format.')
    parser.add_argument('--sort-bam', dest='sort_bam', action="store_true")

    parser.add_argument('--output-unmapped-reads', dest='output_unmapped_reads',
                        help='Additional output file with unmapped reads (single-end).')
    parser.add_argument('--output-unmapped-reads-l', dest='output_unmapped_reads_l',
                        help='File name for unmapped reads (left, paired-end).')
    parser.add_argument('--output-unmapped-reads-r', dest='output_unmapped_reads_r',
                        help='File name for unmapped reads (right, paired-end).')

    parser.add_argument('--output-suppressed-reads', dest='output_suppressed_reads',
                        help='Additional output file with suppressed reads (single-end).')
    parser.add_argument('--output-suppressed-reads-l', dest='output_suppressed_reads_l',
                        help='File name for suppressed reads (left, paired-end).')
    parser.add_argument('--output-suppressed-reads-r', dest='output_suppressed_reads_r',
                        help='File name for suppressed reads (right, paired-end).')
    parser.add_argument('--stdout', dest='output_stdout',
                        help='File name for the standard output of bismark.')

    parser.add_argument('--single-paired', dest='single_paired',
                        help='The single-end reads file in Sanger FASTQ or FASTA format.')

    parser.add_argument('--fastq', action='store_true', help='Query filetype is in FASTQ format')
    parser.add_argument('--fasta', action='store_true', help='Query filetype is in FASTA format')
    parser.add_argument('--phred64-quals', dest='phred64', action="store_true")
    parser.add_argument('--non-directional', dest='non_directional', action="store_true")
    parser.add_argument('--pbat', dest='pbat', action="store_true")

    parser.add_argument('--skip-reads', dest='skip_reads', type=int)
    parser.add_argument('--score-min', dest='score_min', type=str)
    parser.add_argument('--qupto', type=int)

    # paired end options
    parser.add_argument('-I', '--minins', dest='min_insert')
    parser.add_argument('-X', '--maxins', dest='max_insert')
    parser.add_argument('--no-mixed', dest='no_mixed', action="store_true")
    parser.add_argument('--no-discordant', dest='no_discordant', action="store_true")

    # parse general options
    # default 20
    parser.add_argument('--seed-len', dest='seed_len', type=int)
    # default 15
    parser.add_argument('--seed-extention-attempts', dest='seed_extention_attempts', type=int)
    # default 0
    parser.add_argument('--seed-mismatches', dest='seed_mismatches', type=int)
    # default 2
    parser.add_argument('--max-reseed', dest='max_reseed', type=int)


    """
    The number of megabytes of memory a given thread is given to store path
    descriptors in --best mode. Best-first search must keep track of many paths
    at once to ensure it is always extending the path with the lowest cumulative
    cost. Bowtie tries to minimize the memory impact of the descriptors, but
    they can still grow very large in some cases. If you receive an error message
    saying that chunk memory has been exhausted in --best mode, try adjusting
    this parameter up to dedicate more memory to the descriptors. Default: 512.
    """
    parser.add_argument('--chunkmbs', type=int, default=512)

    args = parser.parse_args()

    logger = logging.getLogger('bismark_wrapper')
    logger.setLevel(logging.DEBUG)
    ch = logging.StreamHandler(sys.stdout)
    if args.output_stdout:
        ch.setLevel(logging.WARNING)
        handler = logging.FileHandler(args.output_stdout)
        handler.setLevel(logging.DEBUG)
        logger.addHandler(handler)
    else:
        ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)

    # Create bismark index if necessary.
    index_dir = ""
    tmp_index_dir = None
    if args.own_file:
        logger.info("Create a temporary index with the offered files from the user. "
                    "Utilizing the script: bismark_genome_preparation")
        tmp_index_dir = tempfile.mkdtemp()
        index_path = os.path.join(tmp_index_dir, '.'.join(os.path.split(args.own_file)[1].split('.')[:-1]))
        try:
            # Create a hard link pointing to args.own_file named 'index_path'.fa.
            os.symlink(args.own_file, index_path + '.fa')
        except Exception as e:
            if os.path.exists(tmp_index_dir):
                shutil.rmtree(tmp_index_dir)
            stop_err(logger, 'Error in linking the reference database!\n%s' % e)
        # bismark_genome_preparation needs the complete path to the folder in which the database is stored
        cmd_index = ['bismark_genome_preparation', '--bowtie2', tmp_index_dir]
        try:
            logger.info("Generating index with: '%s'", " ".join(cmd_index))
            process = subprocess.Popen(cmd_index, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with process.stdout:
                log_subprocess_output(logger, process.stdout)
            exitcode = process.wait()
            if exitcode != 0:
                raise Exception(process.stderr)
        except Exception as e:
            if os.path.exists(tmp_index_dir):
                shutil.rmtree(tmp_index_dir)
            stop_err(logger, 'Error indexing reference sequence!\n%s' % e)
        index_dir = tmp_index_dir
    else:
        # bowtie path is the path to the index directory and the first path of the index file name
        index_dir = os.path.dirname(args.index_path)

    # Build bismark command
    tmp_bismark_dir = tempfile.mkdtemp()
    output_dir = os.path.join(tmp_bismark_dir, 'results')
    cmd = ['bismark', '--bam', '--gzip', '--temp_dir', tmp_bismark_dir, '-o', output_dir, '--quiet']

    if args.fasta:
        # the query input files (specified as mate1,mate2 or singles) are FastA
        cmd.append('--fasta')
    elif args.fastq:
        cmd.append('--fastq')

    # alignment options
    if args.num_threads > 2:
        # divide num_threads by 2 here since bismark will spawn 2 jobs with -p threads each
        cmd.extend(['-p', str(math.ceil(args.num_threads / 2))])
    if args.seed_mismatches:
        cmd.extend(['-N', str(args.seed_mismatches)])
    if args.seed_len:
        cmd.extend(['-L', str(args.seed_len)])
    if args.seed_extention_attempts:
        cmd.extend(['-D', str(args.seed_extention_attempts)])
    if args.max_reseed:
        cmd.extend(['-R', str(args.max_reseed)])
    if args.no_discordant:
        cmd.append('--no-discordant')
    if args.no_mixed:
        cmd.append('--no-mixed')
    if args.skip_reads:
        cmd.extend(['--skip', str(args.skip_reads)])
    if args.score_min:
        cmd.extend(['--score_min', str(args.score_min)])
    if args.qupto:
        cmd.extend(['--upto', 'args.qupto'])
    if args.phred64:
        cmd.append('--phred64-quals')
    if args.non_directional:
       cmd.append('--non-directional')
    if args.pbat:
       cmd.append('--pbat')
    if args.suppress_header:
        cmd.append('--sam-no-hd')
    if args.output_unmapped_reads or (args.output_unmapped_reads_l and args.output_unmapped_reads_r):
        cmd.append('--un')
    if args.output_suppressed_reads or (args.output_suppressed_reads_l and args.output_suppressed_reads_r):
        cmd.append('--ambiguous')

    cmd.append(index_dir)
    # Set up the reads
    if args.mate_paired:
        # paired-end reads library
        cmd.extend(['-1', args.mate1, '-2', args.mate2, '-I', str(args.min_insert), '-X', str(args.max_insert)])
    else:
        # single paired reads library
        cmd.append(args.single_paired)

    # Run
    try:
        logger.info("Running bismark with: '%s'", " ".join(cmd))
        process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        with process.stdout:
            log_subprocess_output(logger, process.stdout)
        exitcode = process.wait()
        if exitcode != 0:
            raise Exception(process.stderr)
    except Exception as e:
        stop_err(logger, 'Error in running bismark!\n%s' % e)

    # collect and copy output files
    if args.output_report_file:
        output_report_file = open(args.output_report_file, 'w+')
        for line in fileinput.input(glob(os.path.join(output_dir, '*report.txt'))):
            output_report_file.write(line)
        output_report_file.close()

    if args.output_suppressed_reads:
        if glob(os.path.join(output_dir, '*ambiguous_reads.*')):
            shutil.move(glob(os.path.join(output_dir, '*ambiguous_reads.*'))[0], args.output_suppressed_reads)
    if args.output_suppressed_reads_l:
        if glob(os.path.join(output_dir, '*ambiguous_reads_1.*')):
            shutil.move(glob(os.path.join(output_dir, '*ambiguous_reads_1.*'))[0], args.output_suppressed_reads_l)
    if args.output_suppressed_reads_r:
        if glob(os.path.join(output_dir, '*ambiguous_reads_2.*')):
            shutil.move(glob(os.path.join(output_dir, '*ambiguous_reads_2.*'))[0], args.output_suppressed_reads_r)

    if args.output_unmapped_reads:
        if glob(os.path.join(output_dir, '*unmapped_reads.*')):
            shutil.move(glob(os.path.join(output_dir, '*unmapped_reads.*'))[0], args.output_unmapped_reads)
    if args.output_unmapped_reads_l:
        if glob(os.path.join(output_dir, '*unmapped_reads_1.*')):
            shutil.move(glob(os.path.join(output_dir, '*unmapped_reads_1.*'))[0], args.output_unmapped_reads_l)
    if args.output_unmapped_reads_r:
        if glob(os.path.join(output_dir, '*unmapped_reads_2.*')):
            shutil.move(glob(os.path.join(output_dir, '*unmapped_reads_2.*'))[0], args.output_unmapped_reads_r)

    try:
        # merge all bam files
        tmp_res = tempfile.NamedTemporaryFile(dir=output_dir).name

        bam_files = glob(os.path.join(output_dir, '*.bam'))
        if len(bam_files) > 1:
            cmd = ['samtools', 'merge', '-@', str(args.num_threads), '-f', tmp_res, ' '.join(bam_files)]
            logger.info("Merging bams with: '%s'", cmd)
            process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            with process.stdout:
                log_subprocess_output(logger, process.stdout)
            exitcode = process.wait()
            if exitcode != 0:
                raise Exception(process.stderr)
        else:
            tmp_res = bam_files[0]

        bam_path = "%s" % tmp_res

        if os.path.exists(bam_path):
            if args.sort_bam:
                cmd = ['samtools', 'sort', '-@', str(args.num_threads), bam_path, 'sorted_bam']
                logger.info("Sorting bam with: '%s'", cmd)
                process = subprocess.Popen(args=cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
                with process.stdout:
                    log_subprocess_output(logger, process.stdout)
                exitcode = process.wait()
                if exitcode != 0:
                    raise Exception(process.stderr)
                shutil.move('sorted_bam.bam', args.output)
            else:
                shutil.move(bam_path, args.output)
        else:
            stop_err(logger, 'BAM file no found:\n%s' % bam_path)

    except Exception as e:
        stop_err(logger, 'Error in merging bam files!\n%s' % e)

    # Clean up temp dirs
    if tmp_index_dir and os.path.exists(tmp_index_dir):
        shutil.rmtree(tmp_index_dir)
    if os.path.exists(tmp_bismark_dir):
        shutil.rmtree(tmp_bismark_dir)
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)


if __name__ == "__main__": __main__()
