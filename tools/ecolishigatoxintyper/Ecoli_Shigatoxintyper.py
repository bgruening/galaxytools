#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
############################################################################
# Istituto Superiore di Sanita'
# European Union Reference Laboratory (EU-RL) for Escherichia coli, including Verotoxigenic E. coli (VTEC)
# Developer: Arnold Knijn arnold.knijn@iss.it
############################################################################
"""

import argparse
import datetime
import os
import shutil
import subprocess
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), 'scripts'))
try:
    import HTML
except ImportError:
    import html as HTML

TOOL_DIR = os.path.dirname(os.path.abspath(__file__))


def insertFile(filename, report):
    with open(filename) as html_in:
        for line in html_in:
            report.write(line)


def insertFileAsTable(filename, report, hasheader=False, tabclass="table table-rep"):
    with open(filename) as table_in:
        table_data = [[str(col) for col in row.split('\t')] for row in table_in]
    insertTable(table_data, report, hasheader, tabclass)


def insertTable(table_data, report, hasheader=False, tabclass="table table-rep"):
    if hasheader:
        htmlcode = HTML.table(table_data[1:], attribs={'class': tabclass}, header_row=table_data[0])
    else:
        htmlcode = HTML.table(table_data, attribs={'class': tabclass})
    report.write(htmlcode)


def openFileAsTable(filename):
    with open(filename) as table_in:
        table_data = [[str(col).rstrip() for col in row.split('\t')] for row in table_in]
    return table_data


def __main__():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_type', dest='input_type', help='type of input files (FASTQ or Contigs)')
    parser.add_argument('-1', '--input1', dest='input1', help='forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('--input1_ext', dest='input1_ext', help='extension of forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('--input1_name', dest='input1_name', help='name of forward or single-end reads file in Sanger FASTQ format')
    parser.add_argument('-2', '--input2', dest='input2', help='reverse reads file in Sanger FASTQ format')
    parser.add_argument('--input2_ext', dest='input2_ext', help='extension of reverse reads file in Sanger FASTQ format')
    parser.add_argument('--input2_name', dest='input2_name', help='name of reverse reads file in Sanger FASTQ format')
    parser.add_argument('--html1', dest='html1', help='html FASTQC file')
    parser.add_argument('--html1_id', dest='html1_id', help='html FASTQC file id')
    parser.add_argument('--html1_path', dest='html1_path', help='html FASTQC file path')
    parser.add_argument('--html2', dest='html2', help='html FASTQC file')
    parser.add_argument('--html2_id', dest='html2_id', help='html FASTQC file id')
    parser.add_argument('--html2_path', dest='html2_path', help='html FASTQC file path')
    parser.add_argument('--log', dest='logfile', help='log file')
    parser.add_argument('--shigatoxin', dest='shigatoxin', help='Blastn for Shiga toxin')
    parser.add_argument('--output', dest='output', help='output report html file')
    args = parser.parse_args()

    try:
        trimmomatic_jar = subprocess.check_output(["which", "trimmomatic.jar"], text=True).strip()
        if trimmomatic_jar:
            os.symlink(trimmomatic_jar, "trimmomatic.jar")
    except (subprocess.CalledProcessError, FileNotFoundError):
        try:
            trimmomatic_path = subprocess.check_output(["which", "trimmomatic"], text=True).strip()
            trimmomatic_jar_path = os.path.realpath(trimmomatic_path) + ".jar"
            if os.path.exists(trimmomatic_jar_path):
                os.symlink(trimmomatic_jar_path, "trimmomatic.jar")
        except (subprocess.CalledProcessError, FileNotFoundError, OSError):
            pass
    log = open(args.logfile, 'w')
    log.write("E coli Shiga toxin typer v2.0\n\nTool versions\n=============\n")
    if args.input2:
        # PAIRED-END READS, FASTQC
        subprocess.call("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input1 + " -d " + args.html1_path + " -o " + args.html1 + " -t text1 -f " + args.input1_ext + " -j " + args.input1_name + " -e " + "fastqc", shell=True)
        shutil.rmtree(args.html1_path)
        subprocess.call("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input2 + " -d " + args.html2_path + " -o " + args.html2 + " -t text2 -f " + args.input2_ext + " -j " + args.input2_name + " -e " + "fastqc", shell=True)
        shutil.rmtree(args.html2_path)
        log.write(os.popen("fastqc -v").read())
        # TRIMMING
        subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar PE -threads ${GALAXY_SLOTS:-6} -phred33 " + args.input1 + " " + args.input2 + " trimmed1 trimmed1unpaired trimmed2 trimmed2unpaired SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36", shell=True)
        log.write("\nTrimmomatic v0.39\n")
        log.write("parameters: SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:36\n")
        # ASSEMBLY
        subprocess.call(TOOL_DIR + "/scripts/stx_subtype_pe.sh trimmed1 trimmed2 " + TOOL_DIR, shell=True)
        # LOGGING
        log.write("\n\nFiltering\n=========\nduk v20110303\n")
        log.write("parameters: k=23\n")
        log.write("\n\nAssembly\n========\n")
        log.write(os.popen("spades.py -v").read())
        log.write("parameters: --isolate, pe1-ff, pe1-1, pe1-2\n\n")
        log.write("SKESA v.2.3.0\n")
    else:
        if args.input_type == "fastq":
            # SINGLE-END READS, FASTQC
            subprocess.call("python " + TOOL_DIR + "/scripts/rgFastQC.py -i " + args.input1 + " -d " + args.html1_path + " -o " + args.html1 + " -t text1 -f " + args.input1_ext + " -j " + args.input1_name + " -e " + "fastqc", shell=True)
            shutil.rmtree(args.html1_path)
            log.write(os.popen("fastqc -v").read())
            # TRIMMING
            subprocess.call("java ${_JAVA_OPTIONS:--Xmx8G} -jar trimmomatic.jar SE -threads ${GALAXY_SLOTS:-6} -phred33 " + args.input1 + " trimmed1 SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55", shell=True)
            log.write("\nTrimmomatic v0.39\n")
            log.write("parameters: SLIDINGWINDOW:5:20 LEADING:3 TRAILING:3 MINLEN:55\n")
            # ASSEMBLY
            subprocess.call(TOOL_DIR + "/scripts/stx_subtype_se.sh trimmed1 " + TOOL_DIR, shell=True)
            # LOGGING
            log.write("\n\nFiltering\n=========\nduk v20110303\n")
            log.write("parameters: k=23\n")
            log.write("\n\nAssembly\n========\n")
            log.write(os.popen("spades.py -v").read())
            log.write("parameters: careful, iontorrent\n")
            log.write("SKESA v.2.3.0\n")
        else:
            # CONTIGS
            os.symlink(args.input1, "stx.fasta")
    # BLAST
    subprocess.call(TOOL_DIR + "/scripts/stx_subtype_fa.sh stx.fasta " + TOOL_DIR, shell=True)
    # SHIGATOXINTYPER: OUTPUT
    shutil.copy("stx_blastn", args.shigatoxin)
    shigatoxin_typing = openFileAsTable("blastn_shigatoxin_fc")
    # LOGGING
    log.write("\n\nShigatoxinTyper\n===============\n")
    log.write(os.popen("blastn -version").read())
    log.write("parameters: evalue=0.001, strand=both, dust=yes, max_target_seqs=20, perc_identity=95.0\n\n")
    log.write(os.popen("cat " + TOOL_DIR + "/data/ShigatoxinTyping_db.txt").read())
    log.close()
    try:
        report = open(args.output, 'w')
        # write head html
        insertFile(TOOL_DIR + "/scripts/report_head.html", report)
        report.write("<td><h1><em>E coli</em> Shiga toxin typer</h1><h2>Report for %s</h2>%s</td>" % (args.input1_name, datetime.datetime.utcnow().strftime("%Y-%m-%d %H:%M UTC")))
        insertFile(TOOL_DIR + "/scripts/report_head2.html", report)
        # write results
        report.write("<h3>Summary</h3>\n")
        if len(shigatoxin_typing) == 0:
            report.write("No subtype match found")
        else:
            # get corresponding subtypes
            shigatoxin_subtypes = []
            shigatoxin_subtypes_raw = []
            shigatoxin_types = openFileAsTable(TOOL_DIR + "/data/stx_subtypes")
            for subtype in shigatoxin_typing:
                blast_pident_100 = float(subtype[1]) == 100
                if (blast_pident_100):
                    for item in shigatoxin_types:
                        if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                            shigatoxin_subtypes.append(item[1])
                            shigatoxin_subtypes_raw.append(item[1])
            # partial matches
            for subtype in shigatoxin_typing:
                for item in shigatoxin_types:
                    if item[0] == subtype[0] and item[1] not in shigatoxin_subtypes_raw:
                        if item[1][0:4] == "stx1":
                            shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[1])) + ")")
                            shigatoxin_subtypes_raw.append(item[1])
                        if item[1][0:4] == "stx2":
                            shigatoxin_subtypes.append(item[1] + "(" + str(float(subtype[1])) + ")")
                            shigatoxin_subtypes_raw.append(item[1])
            shigatoxin_subtypes.sort()
            str_shigatoxin_subtype = " ".join(shigatoxin_subtypes)
            report.write(str_shigatoxin_subtype)
            if "(" in str_shigatoxin_subtype:
                report.write("<br/><br/>(*) indicates no complete match was found, the number is the percentage of the sequence matched")
        if args.input_type == "fastq":
            report.write("<br/><br/><hr/><h3>Raw data quality check</h3>\n")
            if args.html1_id:
                report.write("<p>FASTQC result forward available in the history</p>\n")
            if args.input2 and args.html2_id:
                report.write("<p>FASTQC result reverse available in the history</p>\n")
        report.write("<br/><hr/><h3>Shiga toxin typing</h3>\n")
        insertFileAsTable("blastn_shigatoxin_fct", report, True)
        # write tail html
        insertFile(TOOL_DIR + "/scripts/report_tail.html", report)
    finally:
        report.close()


if __name__ == "__main__":
    __main__()
