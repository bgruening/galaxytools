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
parser.add_argument( '--operation', dest='operation', help='Operation' )
parser.add_argument( '--expression', dest='expression', default=None, help='Expression' )
parser.add_argument( '--bin_constant', dest='bin_constant', type=int, default=None, help='Constant of type binary integer' )
parser.add_argument( '--float_constant', dest='float_constant', type=float, default=None, help='Constant of type float' )
parser.add_argument( '--jython_script', dest='jython_script', help='Path to the Jython script' )
parser.add_argument( '--output', dest='output', help='Path to the output file' )
parser.add_argument( '--output_datatype', dest='output_datatype', help='Datatype of the output image' )
args = parser.parse_args()

tmp_dir = imagej2_base_utils.get_temp_dir()
# ImageJ expects valid image file extensions, so the Galaxy .dat extension does not
# work for some features.  The following creates a symlink with an appropriate file
# extension that points to the Galaxy dataset.  This symlink is used by ImageJ.
tmp_input_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.input, args.input_datatype )
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
cmd += ' %s' % args.operation
# Handle the expression, which must be enclosed in " if not None.
if args.expression in [ None, 'None' ]:
    cmd += ' None'
else:
    cmd += ' "%s"' % args.expression
cmd += imagej2_base_utils.handle_none_type( args.bin_constant, val_type='int' )
cmd += imagej2_base_utils.handle_none_type( args.float_constant )
cmd += ' %s' % tmp_output_path
cmd += ' %s' % args.output_datatype

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

# Save the output image.
shutil.move( tmp_output_path, args.output )
imagej2_base_utils.cleanup_before_exit( tmp_dir )
