#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import tempfile
import imagej2_base_utils

parser = argparse.ArgumentParser()
parser.add_argument( '--input', dest='input', help='Path to the input file' )
parser.add_argument( '--input_datatype', dest='input_datatype', help='Datatype of the input image' )
parser.add_argument( '--black_background', dest='black_background', help='Black background' )
parser.add_argument( '--size', dest='size', help='Size (pixel^2)' )
parser.add_argument( '--circularity_min', dest='circularity_min', type=float, help='Circularity minimum' )
parser.add_argument( '--circularity_max', dest='circularity_max', type=float, help='Circularity maximum' )
parser.add_argument( '--show', dest='show', help='Show' )
parser.add_argument( '--display_results', dest='display_results', help='Display results' )
parser.add_argument( '--all_results', dest='all_results', help='All results' )
parser.add_argument( '--exclude_edges', dest='exclude_edges', help='Exclude edges' )
parser.add_argument( '--include_holes', dest='include_holes', help='Include holes' )
parser.add_argument( '--jython_script', dest='jython_script', help='Path to the Jython script' )
parser.add_argument( '--results', dest='results', default=None, help='Path to the output results file' )
parser.add_argument( '--output', dest='output', default=None, help='Path to the output image file' )
parser.add_argument( '--output_datatype', dest='output_datatype', default='data',  help='Datatype of the output image' )
args = parser.parse_args()

tmp_dir = imagej2_base_utils.get_temp_dir()
# ImageJ expects valid image file extensions, so the Galaxy .dat extension does not
# work for some features.  The following creates a symlink with an appropriate file
# extension that points to the Galaxy dataset.  This symlink is used by ImageJ.
tmp_input_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.input, args.input_datatype )
if args.output is None:
    tmp_output_path = None
else:
    tmp_output_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.output_datatype )

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
cmd += ' %s' % args.black_background
cmd += ' %s' % args.size
cmd += ' %.3f' % args.circularity_min
cmd += ' %.3f' % args.circularity_max
cmd += ' %s' % args.show
cmd += ' %s' % args.display_results
cmd += '%s' % imagej2_base_utils.handle_none_type( args.all_results, val_type='str' )
cmd += ' %s' % args.exclude_edges
cmd += ' %s' % args.include_holes
cmd += '%s' % imagej2_base_utils.handle_none_type( tmp_output_path, val_type='str' )
cmd += ' %s' % args.output_datatype
cmd += '%s' % imagej2_base_utils.handle_none_type( args.results, val_type='str' )

# Run the command.
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

if tmp_output_path is not None:
    # Save the output image.
    shutil.move( tmp_output_path, args.output )
imagej2_base_utils.cleanup_before_exit( tmp_dir )
