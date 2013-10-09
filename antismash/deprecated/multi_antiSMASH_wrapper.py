#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys, subprocess, commands
import random, shutil
import zipfile
from Bio import SeqIO
import tempfile

blastdbpath = '/home/galaxy/bin/antismash-1.1.0/db'
pfamdbpath = '/home/galaxy/bin/antismash-1.1.0/db'
antismash_path = '/home/galaxy/bin/antismash-1.1.0/antismash.py'


def anitSMASH(args, sequence_path, glimmer_path):
    #./antismash.py Tue6071_genome.fasta --geneclustertypes 1 --fullhmm y
    print sequence_path
    rint = random.randint(1,10000000)
    tmp_dir = '/tmp/galaxy_%s' % rint
    os.mkdir(tmp_dir)
    os.mkdir(os.path.join( tmp_dir, 'geneprediction' ))
    os.chdir(tmp_dir)
    print 'IFORMAT:',args.iformat
    if args.iformat == 'fasta':
        new_input_path = os.path.join(tmp_dir, os.path.basename(sequence_path) + '.fasta')
    elif args.iformat == 'genbank':
        new_input_path = os.path.join(tmp_dir, os.path.basename(sequence_path) + '.gbk')
    else:
        new_input_path = os.path.join(tmp_dir, os.path.basename(sequence_path) + '.embl')


    # try to generate the same name as in antismash.py
    genomename = ".".join( (os.path.basename(sequence_path) + '.'+ args.iformat).split(".")[:-1] )
    for i in """!"#$%&()*+,./:;=>?@[]^`{|}'""":
        genomename = genomename.replace(i,"")
    result_path = os.path.join( tmp_dir, genomename )

    shutil.copy(sequence_path, new_input_path )

    if args.eukaryotic:
        taxon = '--taxon e'
    else:
        taxon = '--taxon p'

    if args.clusterblast:
        clusterblast = '--clusterblast y'
    else:
        clusterblast = '--clusterblast n'

    if args.smcogs:
        smcogs = '--smcogs y'
    else:
        smcogs = '--smcogs n'

    if args.fullhmm:
        fullhmm = '--fullhmm y'
    else:
        fullhmm = '--fullhmm n'

    if args.fullblast:
        fullblast = '--fullblast y'
    else:
        fullblast = '--fullblast n'

    h = [antismash_path, new_input_path, 
        '--geneclustertypes %s' % args.geneclustertypes, 
        taxon, 
        clusterblast, 
        smcogs, 
        fullhmm,
        fullblast,
        '--blastdbpath %s' % blastdbpath,
        '--pfamdbpath %s' % pfamdbpath,
        '--cores 10',
        ]
    if args.iformat == 'fasta':
        h.append('--glimmer_prediction %s' % glimmer_path)

    a = ' '.join(h)

    subprocess.call(a, shell=True)
    print a
    embl_path = os.path.join(result_path, '%s.final.embl' % genomename)
    if os.path.exists(embl_path) and args.embl_path:
        args.embl_path.write( open(embl_path).read() )

    gff_line = "%s\tantismash\t%s\t%s\t%s\t.\t%s\t.\t%s"

    clustername_mapping = {}
    genecluster_path = os.path.join(result_path, 'clusterblast/geneclusters.txt')
    if os.path.exists(genecluster_path):
        for line in open( genecluster_path):
            token = line.split('\t')
            clustername_mapping[token[2]] = token[3]
        old_cluster_name = False
        parent_cluster_name = ""
        child_lines = ""
        starts = []
        ends = []
        for line in open( os.path.join(result_path, 'clusterblast/geneclusterprots.fasta') ):
            if line.startswith('>'):
                columns = line[1:].strip().split('|')
                args.geneclusterprots.write( line.replace('|%s|' % columns[1], '|%s|%s|' % (columns[1],clustername_mapping[columns[1]])) )

                """
                "Handlampe" kann man zu einer Person sagen wenn es nicht mehr zum Armleuchter reicht. Es ist also ein abgeschwächte Form für "Depp".
                """
                if old_cluster_name == False:
                    # inital clustername setting
                    # try to catch the line in which the name changes
                    old_cluster_name = columns[1]
                    parent_cluster_name = columns[1]

                seqname = columns[0]
                strand = columns[3]
                (start, end) = columns[2].split('-')
                starts.append(int(start))
                ends.append(int(end))

                if old_cluster_name != columns[1]:
                    # begin of a new cluster, caclulate the parent line for the previous clusterproteins
                    # TODO: maybe the assumption is not correct, to take the smallest start and the largest end
                    starts.sort()
                    ends.sort()
                    genecluster_start = starts[0]
                    genecluster_end = ends[-1]
                    group = 'ID=%s;Name=%s;Clustertype=%s;Note=%s\n' % (parent_cluster_name, columns[4], clustername_mapping[columns[1]], columns[5])
                    oline = gff_line % (
                            seqname, 'genecluster', genecluster_start, genecluster_end, strand, group
                    )
                    args.geneclusterprots_gff.write(oline + child_lines)
                    starts = []
                    ends = []
                    old_cluster_name = columns[1]
                    parent_cluster_name = columns[1]
                    child_lines = ''

                group = 'ID=%s;Parent=%s;Name=%s;Clustertype=%s;Note=%s\n' % (columns[4], columns[1], columns[4], clustername_mapping[columns[1]], columns[5])
                child_lines += gff_line % (
                    seqname, 'gene', start, end, strand, group
                )

            else:
                args.geneclusterprots.write( line )

        starts.sort()
        ends.sort()
        genecluster_start = starts[0]
        genecluster_end = ends[-1]
        group = 'ID=%s;Name=%s;Clustertype=%s;Note=%s\n' % (parent_cluster_name, columns[4], clustername_mapping[columns[1]], columns[5])
        oline = gff_line % (
                seqname, 'genecluster', genecluster_start, genecluster_end, strand, group
        )
        args.geneclusterprots_gff.write(oline + child_lines)

    # remove tmp directory
    shutil.rmtree(tmp_dir)


def arg_parse():
    import argparse
    parser = argparse.ArgumentParser(prog = 'antiSMASH-Wrapper')
    parser.add_argument('--version', action='version', version='%(prog)s 0.01')
    parser.add_argument('--geneclustertypes')
    parser.add_argument('--clusterblast', action='store_true')
    parser.add_argument('--eukaryotic', action='store_true')
    parser.add_argument('--fullhmm', action='store_true')
    parser.add_argument('--smcogs', action='store_true')
    parser.add_argument('--fullblast', action='store_true')

    parser.add_argument('--input', '-i', help='FASTA Sequence File')
    parser.add_argument('--glimmer_prediction', help='Glimmer Prediction File')

    parser.add_argument('--input_format', '-f', dest="iformat", help='input format, if the input is a fasta file you need a glimmer file as additional input')

    parser.add_argument('--zip', help='output: all files as zip file')
    parser.add_argument('--html_file', help='output: the path to the index html file')
    parser.add_argument('--html_path', help='output: the path to the output html dir')
    parser.add_argument('--embl_path', help='output: the path to the embl output file', type=argparse.FileType('w'))
    parser.add_argument('--geneclusterprots', help='output: Genecluster Fasta File', type=argparse.FileType('w'))
    parser.add_argument('--geneclusterprots_gff', help='output: Genecluster GFF', type=argparse.FileType('w'))


    args = parser.parse_args()
    return args


def extract_glimmerHMM_results(sequence_id, glimmerHMM):
    """
        Returns the given annoations from a sequence_id as complete GFF3 formated string.
    """
    result = "##gff-version 3\n"
    found_annoations = False
    for line in open(glimmerHMM):
        if line.startswith(sequence_id + '\t'):
            result += line
            found_annoations = True
    if found_annoations:
        return result
    else:
        return False

def extract_glimmer_results(sequence_id, glimmer):
    """
        Returns the given annoations from a sequence_id as glimmer result file.
    """
    result = ''
    found_annoations = False
    for line in open(glimmer):
        # that construct matches the pseudo fasta header and the corresponding orfs
        if line.split()[0].startswith(sequence_id):
            result += line
            found_annoations = True
    if found_annoations:
        return result
    else:
        return False


if __name__ == '__main__':
    args = arg_parse()
    geneclusterprots = ""
    temp_dir = tempfile.mkdtemp()
    args.geneclusterprots_gff.write("##gff-version 3\n")



    for record in SeqIO.parse(open(args.input, "rU"), args.iformat) :
        # create two tem files and fill it with the information from one sequence
        record.id = record.id or record.name
        tmp_sequence = os.path.join(temp_dir, record.id) #fasta -file
        SeqIO.write(record, open(tmp_sequence, 'w+'), args.iformat)

        tmp_glimmer = ''
        if args.iformat == 'fasta':
            # We need to annotate the sequences manually
            if args.eukaryotic:
                #print os.path.join(temp_dir, record.id + '.gff')
                tmp_glimmer = os.path.join(temp_dir, record.id + '.gff')
                annotations = extract_glimmerHMM_results(record.id, args.glimmer_prediction)
                #print annotations, ':', record.id, args.glimmer_prediction
                if not annotations:
                    continue
                open(tmp_glimmer, 'w+').write( annotations )
            else:
                tmp_glimmer = os.path.join(temp_dir, record.id + '.tsv')
                annotations = extract_glimmer_results(record.id, args.glimmer_prediction)
                if not annotations:
                    continue
                open(tmp_glimmer, 'w+').write( annotations )

        anitSMASH(args, tmp_sequence, tmp_glimmer)
    shutil.rmtree(temp_dir)












