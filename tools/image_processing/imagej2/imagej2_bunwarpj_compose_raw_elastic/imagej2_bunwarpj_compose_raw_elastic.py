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
parser.add_argument( '--source_elastic_transformation', dest='source_elastic_transformation', help='Direct source transformation matrix' )
parser.add_argument( '--target_raw_transformation', dest='target_raw_transformation', help='Inverse target transformation matrix' )
parser.add_argument( '--max_heap_size_type', dest='max_heap_size_type', help='Type (default or megabytes) of max_heap_size value' )
parser.add_argument( '--max_heap_size', dest='max_heap_size', help='Maximum size of the memory allocation pool used by the JVM.' )
parser.add_argument( '--output', dest='output', help='Warping index' )

args = parser.parse_args()

tmp_dir = imagej2_base_utils.get_temp_dir()
source_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.source_image, args.source_image_format )
target_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.target_image, args.target_image_format )
source_elastic_transformation_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.source_elastic_transformation, 'txt' )
target_raw_transformation_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.target_raw_transformation, 'txt' )

# Set the size of the memory allocation pool used by the JVM.
memory_size = imagej2_base_utils.get_max_heap_size_value( args.max_heap_size_type, args.max_heap_size )

# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open( tmp_out, 'wb' )
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open( tmp_err, 'wb' )

# Build the command line to compose the raw and elastic transformations.
cmd = imagej2_base_utils.get_base_cmd_bunwarpj( memory_size )
if cmd is None:
    imagej2_base_utils.stop_err( "bUnwarpJ not found!" )
cmd += ' -compose_raw_elastic'
# Target is sent before source.
cmd += ' %s' % target_image_path
cmd += ' %s' % source_image_path
cmd += ' %s' % target_raw_transformation_path
cmd += ' %s' % source_elastic_transformation_path
cmd += ' %s' % args.output

# Compose the raw and elastic transformations into another raw transformation using bUnwarpJ.
proc = subprocess.Popen( args=cmd, stderr=subprocess.PIPE, stdout=subprocess.PIPE, shell=True )
rc = proc.wait()
if rc != 0:
    error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
    imagej2_base_utils.stop_err( error_message )

imagej2_base_utils.cleanup_before_exit( tmp_dir )
