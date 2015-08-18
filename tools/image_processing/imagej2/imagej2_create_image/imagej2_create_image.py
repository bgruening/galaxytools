#!/usr/bin/env python
import argparse
import shutil
import subprocess
import tempfile
import imagej2_base_utils

if __name__=="__main__":
    # Parse Command Line.
    parser = argparse.ArgumentParser()
    parser.add_argument( '--width', dest='width', type=int, help='Image width in pixels' )
    parser.add_argument( '--height', dest='height', type=int, help='Image height in pixels' )
    parser.add_argument( '--depth', dest='depth', type=int, help='Image depth (specifies the number of stack slices)' )
    parser.add_argument( '--image_type', dest='image_type', help='Image type' )
    parser.add_argument( '--image_title', dest='image_title', default='', help='Image title' )
    parser.add_argument( '--output_datatype', dest='output_datatype', help='Output image format' )
    parser.add_argument( '--jython_script', dest='jython_script', help='Path to the Jython script' )
    parser.add_argument( '--out_fname', help='Path to the output file' )
    args = parser.parse_args()

    tmp_dir = imagej2_base_utils.get_temp_dir()
    tmp_image_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.output_datatype )

    # Define command response buffers.
    tmp_out = tempfile.NamedTemporaryFile().name
    tmp_stdout = open( tmp_out, 'wb' )
    tmp_err = tempfile.NamedTemporaryFile().name
    tmp_stderr = open( tmp_err, 'wb' )
    # Build the command line.
    cmd = imagej2_base_utils.get_base_command_imagej2( None, jython_script=args.jython_script )
    if cmd is None:
        imagej2_base_utils.stop_err( "ImageJ not found!" )
    cmd += ' %s %d %d %d %s %s' % ( args.image_title, args.width, args.height, args.depth, args.image_type, tmp_image_path )
    proc = subprocess.Popen( args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True )
    rc = proc.wait()
    if rc != 0:
        error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
        imagej2_base_utils.stop_err( error_message )
    shutil.move( tmp_image_path, args.out_fname )
    imagej2_base_utils.cleanup_before_exit( tmp_dir )
