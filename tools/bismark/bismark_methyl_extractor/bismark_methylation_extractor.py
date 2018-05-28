#!/usr/bin/env python

import argparse, os, shutil, subprocess, sys, tempfile, fileinput
import zipfile
import re
from glob import glob

def stop_err( msg ):
    sys.stderr.write( "%s\n" % msg )
    sys.exit()

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
        """
            Create a hard link pointing to genome_file named 'genome_path'.fa.
        """
        os.symlink(genome_file, genome_path + '.fa')
    except Exception, e:
        if os.path.exists(tmp_genome_dir):
            shutil.rmtree(tmp_genome_dir)
        stop_err('Error in linking the reference database.\n' + str(e))
    return tmp_genome_dir

def __main__():
    #Parse Command Line
    parser = argparse.ArgumentParser(description='Wrapper for the bismark methylation caller.')

    # input options
    parser.add_argument( '--bismark_path', dest='bismark_path', help='Path to the bismark perl scripts' )

    parser.add_argument( '--infile', help='Input file in SAM or BAM format.' )
    parser.add_argument( '--single-end', dest='single_end', action="store_true" )
    parser.add_argument( '--paired-end', dest='paired_end', action="store_true" )

    parser.add_argument('--splitting_report', dest='splitting_report')
    parser.add_argument('--mbias_report', dest='mbias_report')
    parser.add_argument('--cytosine_report', dest="cytosine_report")
    parser.add_argument('--genome_file', dest="genome_file")
    parser.add_argument('--cx_context', action="store_true" )

    parser.add_argument( '--comprehensive', action="store_true" )
    parser.add_argument( '--merge-non-cpg', dest='merge_non_cpg', action="store_true" )
    parser.add_argument( '--no-overlap', dest='no_overlap', action="store_true" )
    parser.add_argument( '--compress' )
    parser.add_argument('--ignore', dest='ignore', type=int)
    parser.add_argument('--ignore_r2', dest='ignore_r2', type=int)
    parser.add_argument('--ignore_3prime', dest='ignore_3prime', type=int)
    parser.add_argument('--ignore_3prime_r2', dest='ignore_3prime_r2', type=int)

    args = parser.parse_args()

    # Build methylation extractor command
    output_dir = tempfile.mkdtemp()
    cmd = 'bismark_methylation_extractor --no_header -o %s %s %s'
    if args.bismark_path:
        # add the path to the bismark perl scripts, that is needed for galaxy
        cmd = os.path.join(args.bismark_path, cmd)

    # Set up all options
    additional_opts = ''
    if args.single_end:
        additional_opts += ' --single-end '
    else:
        additional_opts += ' --paired-end '
    if args.no_overlap:
        additional_opts += ' --no_overlap '
    if args.ignore:
        additional_opts += ' --ignore %s ' % args.ignore
    if args.ignore_r2:
        additional_opts += ' --ignore_r2 %s ' % args.ignore_r2
    if args.ignore_3prime:
        additional_opts += ' --ignore_3prime %s ' % args.ignore_3prime
    if args.ignore_3prime_r2:
        additional_opts += ' --ignore_3prime_r2 %s ' % args.ignore_3prime_r2
    if args.comprehensive:
        additional_opts += ' --comprehensive '
    if args.merge_non_cpg:
        additional_opts += ' --merge_non_CpG '
    if args.splitting_report:
        additional_opts += ' --report '
    if args.cytosine_report:
        tmp_genome_dir = build_genome_dir(args.genome_file)
        if args.cx_context:
            additional_opts += ' --bedgraph --CX_context --cytosine_report --CX_context --genome_folder %s ' % tmp_genome_dir
        else:
            additional_opts += ' --bedgraph --cytosine_report --genome_folder %s ' % tmp_genome_dir


    #detect BAM file, use samtools view if it is a bam file
    f = open (args.infile, 'rb')
    sig = f.read(4)
    f.close()
    if sig == '\x1f\x8b\x08\x04' :
        #cmd = cmd % (output_dir, additional_opts, '-')
        new_infilename = os.path.join(output_dir, 'submitted_bs_mapped_reads.sam')
        new_sam = open(new_infilename, 'wb')
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open(tmp_err, 'wb')
        proc = subprocess.Popen(['samtools', 'view', args.infile], stdout=new_sam, stderr=tmp_stderr)
        new_sam.close()
        tmp_stderr.close()
        if os.stat(tmp_err).st_size != 0:
            tmp_sterr = open(tmp_err, 'rb')
            error_msg = tmp_sterr.read()
            tmp_sterr.close()
            sys.exit("error: %s" % error_msg)
        cmd = cmd % (output_dir, additional_opts, new_infilename)
    else:
        cmd = cmd % (output_dir, additional_opts, args.infile)

    # Run
    try:
        tmp_out = tempfile.NamedTemporaryFile().name
        tmp_stdout = open( tmp_out, 'wb' )
        tmp_err = tempfile.NamedTemporaryFile().name
        tmp_stderr = open( tmp_err, 'wb' )
        proc = subprocess.Popen( args=cmd, shell=True, cwd=".", stdout=tmp_stdout, stderr=tmp_stderr )
        returncode = proc.wait()
        tmp_stderr.close()
        # get stderr, allowing for case where it's very large
        tmp_stderr = open( tmp_err, 'rb' )
        stderr = ''
        buffsize = 1048576
        try:
            while True:
                stderr += tmp_stderr.read( buffsize )
                if not stderr or len( stderr ) % buffsize != 0:
                    break
        except OverflowError:
            pass
        tmp_stdout.close()
        tmp_stderr.close()
        if returncode != 0:
            raise Exception(stderr)
            
        # TODO: look for errors in program output.
    except Exception, e:
        stop_err( 'Error in bismark methylation extractor:\n' + str( e ) )

    # collect and copy output files
    if args.compress:
        zipper(output_dir, args.compress)

    # cytosine report
    if args.cytosine_report:
        if args.cx_context:
            shutil.move( glob(os.path.join( output_dir, '*CX_report.txt'))[0], args.cytosine_report )
        else:
            shutil.move(glob(os.path.join(output_dir, '*CpG_report.txt'))[0], args.cytosine_report)
    # splitting report
    if args.splitting_report:
        shutil.move( glob(os.path.join( output_dir, '*_splitting_report.txt'))[0], args.splitting_report )
    if args.mbias_report:
        shutil.move(glob(os.path.join(output_dir, '*M-bias.txt'))[0], args.mbias_report)


    #Clean up temp dirs
    if os.path.exists( output_dir ):
         shutil.rmtree( output_dir )
    if args.cytosine_report:
        if os.path.exists( tmp_genome_dir ):
            shutil.rmtree( tmp_genome_dir )

if __name__=="__main__": __main__()
