#!/usr/bin/env python
import argparse
import subprocess
import tempfile
import imagej2_base_utils

# Parse Command Line.
parser = argparse.ArgumentParser()
parser.add_argument( '--source_image', dest='source_image', help='Source image' )
parser.add_argument( '--source_image_format', dest='source_image_format', help='Source image format' )
parser.add_argument( '--target_image', dest='target_image', help='Target image' )
parser.add_argument( '--target_image_format', dest='target_image_format', help='Target image format' )
parser.add_argument( '--input_elastic_transformation', dest='input_elastic_transformation', help='Input elastic transformation matrix' )
parser.add_argument( '--image_size_factor', dest='image_size_factor', type=float, help='Image size factor' )
parser.add_argument( '--output', dest='output', help='Warping index' )

args = parser.parse_args()

tmp_dir = imagej2_base_utils.get_temp_dir()
source_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.source_image, args.source_image_format )
target_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.target_image, args.target_image_format )
input_elastic_transformation_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.input_elastic_transformation, 'txt' )

# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open( tmp_out, 'wb' )
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open( tmp_err, 'wb' )

def is_power2( val ):
    if val < 0:
        return False
    if val < 1:
        val = 1.0 / val
    val = int( val )
    return ( ( val & ( val - 1 ) ) == 0 )

# Build the command line to adapt the transformation.
cmd = imagej2_base_utils.get_base_cmd_bunwarpj( None )
if cmd is None:
    imagej2_base_utils.stop_err( "bUnwarpJ not found!" )
cmd += ' -adapt_transform'

# Make sure the value of image_size_factor is a power of 2 (positive or negative).
if is_power2( args.image_size_factor ):
    image_size_factor = args.image_size_factor
else:
    msg = "Image size factor must be a positive or negative power of 2 (0.25, 0.5, 2, 4, 8, etc)."
    imagej2_base_utils.stop_err( msg )

# Target is sent before source.
cmd += ' %s' % target_image_path
cmd += ' %s' % source_image_path
cmd += ' %s' % input_elastic_transformation_path
cmd += ' %s' % args.output
cmd += ' %2.f' % image_size_factor

# Adapt the transformation based on the image size factor using bUnwarpJ.
proc = subprocess.Popen( args=cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True )
rc = proc.wait()
if rc != 0:
    error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
    imagej2_base_utils.stop_err( error_message )

imagej2_base_utils.cleanup_before_exit( tmp_dir )
