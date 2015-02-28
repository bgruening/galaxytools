#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os, sys, subprocess, commands
import random, shutil
import zipfile


blastdbpath = '/home/galaxy/bin/antismash-1.1.0/db'
pfamdbpath = '/home/galaxy/bin/antismash-1.1.0/db'
antismash_path = '/home/galaxy/bin/antismash-1.1.0/antismash.py'


def zipper(dir, zip_file):
    zip = zipfile.ZipFile(zip_file, 'w', compression=zipfile.ZIP_DEFLATED)
    root_len = len(os.path.abspath(dir))
    for root, dirs, files in os.walk(dir):
        # only inlcude the result directory
        # assumption, each galaxy file and so the result directory starts with dataset_xxx
        if root.find('dataset_') != -1:
            archive_root = os.path.abspath(root)[root_len:]
            for f in files:
                fullpath = os.path.join(root, f)
                archive_name = os.path.join(archive_root, f)
                zip.write(fullpath, archive_name, zipfile.ZIP_DEFLATED)
    zip.close()
    return zip_file


def anitSMASH(args):
    #./antismash.py Tue6071_genome.fasta --geneclustertypes 1 --fullhmm y
    rint = random.randint(1,10000000)
    tmp_dir = '/tmp/galaxy_%s' % rint
    os.mkdir(tmp_dir)
    os.mkdir(os.path.join( tmp_dir, 'geneprediction' ))
    os.chdir(tmp_dir)
    new_input_path = os.path.join(tmp_dir, os.path.basename(args.input) + '.fasta')

    # try to generate the same name as in antismash.py
    genomename = ".".join( (os.path.basename(args.input) + '.fasta').split(".")[:-1] )
    for i in """!"#$%&()*+,./:;=>?@[]^`{|}'""":
        genomename = genomename.replace(i,"")
    result_path = os.path.join( tmp_dir, genomename )

    shutil.copy(args.input, new_input_path )

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
        '--glimmer_prediction %s' % args.glimmer_prediction,
        '--blastdbpath %s' % blastdbpath,
        '--pfamdbpath %s' % pfamdbpath,
        '--cores 10',
        ]
    a = ' '.join(h)
    subprocess.call(a, shell=True)


    shutil.copy(os.path.join(result_path, '%s.final.embl' % genomename), args.embl_path)

    clustername_mapping = {}
    for line in open( os.path.join(result_path, 'clusterblast/geneclusters.txt') ):
        token = line.split('\t')
        clustername_mapping[token[2]] = token[3]

    for line in open( os.path.join(result_path, 'clusterblast/geneclusterprots.fasta') ):
        if line.startswith('>'):
            for k,v in clustername_mapping.items():
                if '|%s|' % k in line:
                    args.geneclusterprots.write( line.replace('|%s|' % k, '|%s|%s|' % (k,v)) )
        else:
            args.geneclusterprots.write( line )

    zipper(result_path, args.zip)

    # html output
    shutil.copy( os.path.join(result_path, 'display.xhtml'), args.html_file)
    os.mkdir( args.html_path )
    html_dest_path = os.path.join(args.html_path, 'html/')
    images_dest_path = os.path.join(args.html_path, 'images/')
    svg_dest_path = os.path.join(args.html_path, 'svg/')
    substrspecs_dest_path = os.path.join(args.html_path, 'substrspecs/')
    shutil.copytree( os.path.join(result_path, 'html/'), html_dest_path)
    shutil.copytree( os.path.join(result_path, 'images/'), images_dest_path)
    shutil.copytree( os.path.join(result_path, 'svg/'), svg_dest_path)
    shutil.copytree( os.path.join(result_path, 'substrspecs/'), substrspecs_dest_path)
    shutil.copy( os.path.join(result_path, 'jquery.svg.js'), args.html_path ) 
    shutil.copy( os.path.join(result_path, 'jquery.svgdom.js'), args.html_path ) 
    shutil.copy( os.path.join(result_path, 'jquery-1.4.2.min.js'), args.html_path ) 
    shutil.copy( os.path.join(result_path, 'style.css'), args.html_path ) 

    # remove tmp directory
    shutil.rmtree(tmp_dir)


def arg_parse():
    import argparse
    parser = argparse.ArgumentParser(prog = 'antiSMASH-Wrapper')
    parser.add_argument('--version', action='version', version='%(prog)s 0.01')
    parser.add_argument('--geneclustertypes',
                   help='Fingerprint Type, currently FP2, FP3, FP4')
    parser.add_argument('--clusterblast', action='store_true')
    parser.add_argument('--eukaryotic', action='store_true')
    parser.add_argument('--fullhmm', action='store_true')
    parser.add_argument('--smcogs', action='store_true')
    parser.add_argument('--fullblast', action='store_true')

    parser.add_argument('--input', '-i', help='FASTA Sequence File')
    parser.add_argument('--glimmer_prediction', help='Glimmer Prediction File')

    parser.add_argument('--zip', help='output: all files as zip file')
    parser.add_argument('--html_file', help='output: the path to the index html file')
    parser.add_argument('--html_path', help='output: the path to the output html dir')
    parser.add_argument('--embl_path', help='output: the path to the embl output file')
    parser.add_argument('--geneclusterprots', help='output: Genecluster Fasta File', type=argparse.FileType('w'))

    args = parser.parse_args()
    return args


if __name__ == '__main__':
    args = arg_parse()
    anitSMASH(args)

