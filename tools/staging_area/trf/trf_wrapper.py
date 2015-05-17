#!/usr/bin/env python

"""
Runs TRF on a sequence file.
For use with trf version 4.04

usage: trf_wrapper.py input_file <required parameters> [options]
     Required parameters:
    match
    mismatch
    indels
    match_probability (%, 10-100)
    indel_probability (%, 10-100)
    minimal score (30-150)
    maximum period size (1-2000)
"""

import optparse, os, shutil, subprocess, sys, re
from os.path import abspath

def remove_html(infilename,outfilename):
    infile = open(infilename, "r")
    outfile = open(outfilename, "w")
    opena = re.compile(r'<A\ HREF.*?>')
    closea = re.compile(r'</A>')
    for line in infile:
        line = opena.sub('', line)
        line = closea.sub('', line)
        outfile.write(line)
    infile.close()
    outfile.close()

def stop_err(msg):
    sys.stderr.write("%s\n" % msg)
    sys.exit()

def __main__():
    #Parse arguments
    usage = "usage: %prog input match mismatch indels match_p indel_p min_score max_period [options]"
    parser = optparse.OptionParser(usage)
    parser.add_option( '-f', '--flanking', dest='flanking', action='store_true', help='Flanking sequence' )
    parser.add_option( '-m', '--masked', dest='masked', action='store_true', help='Masked sequence' )
    parser.add_option( '-r', '--redundancy', dest='redundancy', action='store_true', help='No redundance elimination' )
    parser.add_option( '-o', '--datoutput', dest='datoutput', action='store', help='Output data file name' )
    parser.add_option( '-k', '--maskoutput', dest='maskoutput', action='store', help='Output mask file name' )
    parser.add_option( '-t', '--report', dest='report', action='store', help='Report file name' )
    parser.add_option( '-i', '--indices', dest='indices', action='store', help='Indices file name' )
    (opts, arguments) = parser.parse_args()

    # Arguments and options
    if len(arguments) != 8:
        print(usage)
        stop_err('Wrong number of arguments passed')
    output_dat = "%s.%s.%s.%s.%s.%s.%s.%s.dat" % (tuple(arguments))
    output_mask = "%s.%s.%s.%s.%s.%s.%s.%s.mask" % (tuple(arguments))
    output_report = "%s.%s.%s.%s.%s.%s.%s.%s.1.html" % (tuple(arguments))
    output_indices = "%s.%s.%s.%s.%s.%s.%s.%s.1.txt.html" % (tuple(arguments))

    if opts.masked and opts.masked == True:
        arguments.append('-m')
    if opts.flanking and opts.flanking == True:
        arguments.append('-f')
    if opts.redundancy and opts.redundancy  == True:
        arguments.append('-r')
    if opts.datoutput and opts.datoutput  != '':
        output_dat_filename = opts.datoutput
    if opts.maskoutput and opts.maskoutput  != '':
        output_mask_filename = opts.maskoutput
    if opts.report and opts.report  != '':
        output_report_filename = opts.report
    if opts.indices and opts.indices  != '':
        output_indices_filename = opts.indices

    # Run
    cmd = arguments
    cmd.insert(0,"trf")
    # Produce both html and dat files
    cmd.append("-d")
    try:
        proc = subprocess.Popen(args=cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except Exception, err:
        sys.stderr.write("Error invoking command: \n%s\n\n%s\n" % (cmd, err))
        sys.exit(1)
    stdout, stderr = proc.communicate()
    return_code = proc.returncode
    if return_code != 1:
        sys.stdout.write(stdout)
        sys.stderr.write(stderr)
        sys.stderr.write("Return error code %i from command:\n" % return_code)
        sys.stderr.write("%s\n" % cmd)
    else:
        sys.stdout.write(stdout)
        sys.stdout.write(stderr)
    cdir = os.getcwd()
    file_list = os.listdir(cdir)
    try:
        report_file = os.path.basename(output_report)
        remove_html(report_file, output_report_filename)
        shutil.copyfile(os.path.basename(output_dat), output_dat_filename)
#        shutil.copyfile(os.path.basename(output_report), output_report_filename)
        shutil.copyfile(os.path.basename(output_indices), output_indices_filename)
        if opts.masked and opts.masked == True:
            shutil.copyfile(os.path.basename(output_mask), output_mask_filename)
    except Exception, err:
        sys.stderr.write("Error copying output files: \n%s\n" % err)
if __name__=="__main__": __main__()
