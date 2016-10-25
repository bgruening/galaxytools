#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import tempfile
import imagej2_base_utils

if __name__=="__main__":
    # Parse Command Line.
    parser = argparse.ArgumentParser()
    parser.add_argument( '--input', dest='input', help='Path to the input file' )
    parser.add_argument( '--input_datatype', dest='input_datatype', help='Datatype of the input image' )
    parser.add_argument( '--noise', dest='noise', help='Specified noise to add to or remove from the image' )
    parser.add_argument( '--standard_deviation', dest='standard_deviation', type=float, default=None, help='Standard deviation' )
    parser.add_argument( '--radius', dest='radius', type=float, default=None, help='Radius' )
    parser.add_argument( '--threshold', dest='threshold', type=float, default=None, help='Threshold' )
    parser.add_argument( '--which_outliers', dest='which_outliers', default=None, help='Which outliers' )
    parser.add_argument( '--randomj', dest='randomj', default=None, help='RandomJ' )
    parser.add_argument( '--trials', dest='trials', type=float, default=None, help='Trials' )
    parser.add_argument( '--probability', dest='probability', type=float, default=None, help='Probability' )
    parser.add_argument( '--lammbda', dest='lammbda', type=float, default=None, help='Lambda' )
    parser.add_argument( '--order', dest='order', type=int, default=None, help='Order' )
    parser.add_argument( '--mean', dest='mean', type=float, default=None, help='Mean' )
    parser.add_argument( '--sigma', dest='sigma', type=float, default=None, help='Sigma' )
    parser.add_argument( '--min', dest='min', type=float, default=None, help='Min' )
    parser.add_argument( '--max', dest='max', type=float, default=None, help='Max' )
    parser.add_argument( '--insertion', dest='insertion', default=None, help='Insertion' )
    parser.add_argument( '--jython_script', dest='jython_script', help='Path to the Jython script' )
    parser.add_argument( '--output', dest='output', help='Path to the output file' )
    args = parser.parse_args()

    tmp_dir = imagej2_base_utils.get_temp_dir()
    # ImageJ expects valid image file extensions, so the Galaxy .dat extension does not
    # work for some features.  The following creates a symlink with an appropriate file
    # extension that points to the Galaxy dataset.  This symlink is used by ImageJ.
    tmp_input_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.input, args.input_datatype )
    tmp_output_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.input_datatype )

    # Define command response buffers.
    tmp_out = tempfile.NamedTemporaryFile().name
    tmp_stdout = open( tmp_out, 'wb' )
    tmp_err = tempfile.NamedTemporaryFile().name
    tmp_stderr = open( tmp_err, 'wb' )
    # Java writes a lot of stuff to stderr, so we'll specify a file for handling actual errors.
    error_log = tempfile.NamedTemporaryFile( delete=False ).name
    # Build the command line.
    cmd = imagej2_base_utils.get_base_command_imagej2( None, jython_script=args.jython_script )
    if cmd is None:
        imagej2_base_utils.stop_err( "ImageJ not found!" )
    cmd += ' %s' % error_log
    cmd += ' %s' % tmp_input_path
    cmd += ' %s' % args.input_datatype
    cmd += ' %s ' % args.noise
    cmd += imagej2_base_utils.handle_none_type( args.standard_deviation )
    cmd += imagej2_base_utils.handle_none_type( args.radius )
    cmd += imagej2_base_utils.handle_none_type( args.threshold )
    cmd += ' %s' % args.which_outliers
    cmd += ' %s' % args.randomj
    cmd += imagej2_base_utils.handle_none_type( args.trials )
    cmd += imagej2_base_utils.handle_none_type( args.probability )
    cmd += imagej2_base_utils.handle_none_type( args.lammbda )
    cmd += imagej2_base_utils.handle_none_type( args.order, val_type='int' )
    cmd += imagej2_base_utils.handle_none_type( args.mean )
    cmd += imagej2_base_utils.handle_none_type( args.sigma )
    cmd += imagej2_base_utils.handle_none_type( args.min )
    cmd += imagej2_base_utils.handle_none_type( args.max )
    cmd += ' %s' % args.insertion
    cmd += ' %s' % tmp_output_path

    proc = subprocess.Popen( args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True )
    rc = proc.wait()

    # Handle execution errors.
    if rc != 0:
        error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
        imagej2_base_utils.stop_err( error_message )
    # Handle processing errors.
    if os.path.getsize( error_log ) > 0:
        error_message = open( error_log, 'r' ).read()
        imagej2_base_utils.stop_err( error_message )
    # Save the output image.
    shutil.move( tmp_output_path, args.output )
    imagej2_base_utils.cleanup_before_exit( tmp_dir )
