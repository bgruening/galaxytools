#!/usr/bin/env python
import argparse
import shutil
import subprocess
import tempfile
import imagej2_base_utils

# Parse Command Line.
parser = argparse.ArgumentParser()
parser.add_argument( '--source_image', dest='source_image', help='Source image' )
parser.add_argument( '--source_image_format', dest='source_image_format', help='Source image format' )
parser.add_argument( '--target_image', dest='target_image', help='Target image' )
parser.add_argument( '--target_image_format', dest='target_image_format', help='Target image format' )
parser.add_argument( '--elastic_transformation', dest='elastic_transformation', help='Elastic transformation as saved by bUnwarpJ in elastic format' )
parser.add_argument( '--source_out', help='Output source image' )
parser.add_argument( '--source_out_datatype', help='Output registered source image format' )
parser.add_argument( '--jython_script', dest='jython_script', help='Path to the Jython script' )

args = parser.parse_args()

tmp_dir = imagej2_base_utils.get_temp_dir()
source_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.source_image, args.source_image_format )
tmp_source_out_tiff_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, 'tiff' )
tmp_source_out_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.source_out_datatype )
target_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.target_image, args.target_image_format )
elastic_transformation_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.elastic_transformation, 'txt' )

# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open( tmp_out, 'wb' )
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open( tmp_err, 'wb' )

# Build the command line to apply the transformation.
cmd = imagej2_base_utils.get_base_cmd_bunwarpj( None )
if cmd is None:
    imagej2_base_utils.stop_err( "bUnwarpJ not found!" )
cmd += ' -elastic_transform'
# Target is sent before source.
cmd += ' %s' % target_image_path
cmd += ' %s' % source_image_path
cmd += ' %s' % elastic_transformation_path
cmd += ' %s' % tmp_source_out_tiff_path

# Apply the elastic transformation using bUnwarpJ.
proc = subprocess.Popen( args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True )
rc = proc.wait()
if rc != 0:
    error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
    imagej2_base_utils.stop_err( error_message )

# Convert the registered image to the specified output format.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open( tmp_out, 'wb' )
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open( tmp_err, 'wb' )

cmd = imagej2_base_utils.get_base_command_imagej2( None, jython_script=args.jython_script )
if cmd is None:
    imagej2_base_utils.stop_err( "ImageJ not found!" )
cmd += ' %s %s %s' % ( tmp_source_out_tiff_path,
                       args.source_out_datatype,
                       tmp_source_out_path )

proc = subprocess.Popen( args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True )
rc = proc.wait()
if rc != 0:
    error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
    imagej2_base_utils.stop_err( error_message )

# Save the Registered Source Image to the defined output.
shutil.move( tmp_source_out_path, args.source_out )
imagej2_base_utils.cleanup_before_exit( tmp_dir )
