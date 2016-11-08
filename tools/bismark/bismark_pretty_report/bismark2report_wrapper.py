#!/usr/bin/python

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import logging

def cleanup_before_exit(tmp_dir):
    if tmp_dir and os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)

def get_arg():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tool_dir', dest='tool_dir', action='store', nargs=1, metavar='tool_dir', type=str)
    parser.add_argument('--alignment_report', dest='alignment_report', action='store', nargs=1, metavar='alignment_report', type=str)
    parser.add_argument('--dedup_report', dest='dedup_report', action='store', nargs=1, metavar='dedup_report', type=str)
    parser.add_argument('--splitting_report', dest='splitting_report', action='store', nargs=1, metavar='splitting_report', type=str)
    parser.add_argument('--mbias_report', dest='mbias_report', action='store', nargs=1, metavar='mbias_report', type=str)
    parser.add_argument('--nucleotide_report', dest='nucleotide_report', action='store', nargs=1, metavar='nucleotide_report', type=str)
    parser.add_argument('--output_html_report', dest='output_html_report', action='store', nargs=1, metavar='output_html_report', type=str)
    parser.add_argument('--output_html_report_link', dest='output_html_report_link', action='store', nargs=1, metavar='output_html_report_link', type=str)
    parser.add_argument('--log_report', dest='log_report', action='store', nargs=1, metavar='log_report', type=str)
    parser.add_argument('--output_dir', dest='job_dir', action='store', nargs=1, metavar='job_dir', type=str)
    args = parser.parse_args()
    return args

def create_and_write_html_link(job_dir, output_html_report_link, tmp_dir):
    """
    Web browsers don't allow to open a link pointing to the absolute path of a local html file FROM a website page;
    The only way to make such link functional is to integrate the local file inside the web structure of the site.
    Galaxy has been designed such that the child_dir <dataset_[0-9]+_files> of the output_dir is considered as the root
    of the html base tag (i.e <base href="/" /> for the current job running.
    The function proceeds the following steps:
    #1. Extracts the galaxy dir where the output files are stored
    #2. Creating a child dir <dataset_[0-9]+_files> in this output_dir is needed because it is considered as the root of the html base tag
    #   We can extract the exact name of this child dir from the jobs_directory name
    #3. Moves the html file in this child_dir
    """
    output_path_list = output_html_report_link.split('/')
    output_path = '/'.join(output_path_list[0:-1])
    html_root = job_dir.split('/')[-1]
    final_dir = os.path.join(output_path, html_root)
    os.makedirs(final_dir)
    shutil.move(os.path.join(tmp_dir, 'html_report'), os.path.join(final_dir, 'html_report.html'))

    html_report = open(output_html_report_link, 'wb')
    html_report.write('<!DOCTYPE html>\n')
    html_report.write('<head>\n')
    html_report.write('\t<meta http-equiv="content-type" content="text/html; charset=UTF-8">\n')
    html_report.write('\t\t<base href="/" />\n')
    html_report.write('\t\t<a href="html_report.html/" target="_blank">Link to Bismark Pretty Report Page</a>\n')
    html_report.write('</head>')
    html_report.close()

def __main__():
    args = get_arg()

    tmp_dir = tempfile.mkdtemp(prefix='tmp', suffix='')

    if args.log_report:
        logging.basicConfig(level=logging.INFO, filename=args.log_report[0], filemode="a+", format='%(message)s')
    else:
        logging.basicConfig(level=logging.INFO, filename=os.path.join(tmp_dir, 'log_report.txt'), filemode="a+", format='%(message)s')

    alignment_option = '--alignment_report'
    alignment_report = args.alignment_report[0]
    if args.dedup_report:
        dedup_option = '--dedup_report'
        dedup_report = args.dedup_report[0]
    else:
        dedup_option = ''
        dedup_report = ''
    if args.splitting_report:
        splitting_option = '--splitting_report'
        splitting_report = args.splitting_report[0]
    else:
        splitting_option = ''
        splitting_report = ''
    if args.mbias_report:
        mbias_option = '--mbias_report'
        mbias_report = args.mbias_report[0]
    else:
        mbias_option = ''
        mbias_report = ''
    if args.nucleotide_report:
        nucleotide_option = '--nucleotide_report'
        nucleotide_report = args.nucleotide_report[0]
    else:
        nucleotide_option = ''
        nucleotide_report = ''

    proc = subprocess.Popen(['perl', os.path.join(args.tool_dir[0], 'bismark2report'), alignment_option, alignment_report, dedup_option, dedup_report,\
                             splitting_option, splitting_report, mbias_option, mbias_report, nucleotide_option, nucleotide_report,\
                             '--dir', tmp_dir, '--output', 'html_report'],\
                             stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    proc_out, proc_err = proc.communicate()

    cmd = 'perl %s %s %s %s %s %s %s %s %s %s %s --output html_report --dir %s'\
          % (os.path.join(args.tool_dir[0], 'bismark2report'), alignment_option, alignment_report, dedup_option, dedup_report,\
             splitting_option, splitting_report, mbias_option, mbias_report, nucleotide_option, nucleotide_report, tmp_dir)

    logging.info('COMMAND LINE:\n\n%s' % cmd)
    logging.info("__________________________________________________________________\n")
    logging.info("BISMARK PRETTY REPORT STDOUT:\n\n%s" % proc_out)
    if proc_err:
        logging.critical("__________________________________________________________________\n")
        logging.critical("BISMARK PRETTY REPORT ERROR:\n\n%s" % proc_err)
        sys.exit("Bismark pretty report crashed with the folowing error message:\n%s" % proc_err)

    if args.output_html_report:
        shutil.copy(os.path.join(tmp_dir, 'html_report'), args.output_html_report[0])

    #This function writes a link towards the Bismark html page inside an html file.
    #This is needed because the direct visualization of the Bismark html report via Galaxy is ugly
    create_and_write_html_link(args.job_dir[0], args.output_html_report_link[0], tmp_dir)

    cleanup_before_exit(tmp_dir)

if __name__=="__main__": __main__()