#!/usr/bin/env python
import argparse
import os
import shutil
import subprocess
import tempfile
import imagej2_base_utils

# Parse Command Line.
parser = argparse.ArgumentParser()
parser.add_argument( '--source_image', dest='source_image', help='Source image' )
parser.add_argument( '--source_image_format', dest='source_image_format', help='Source image format' )
parser.add_argument( '--source_mask', dest='source_mask', default=None, help='Source mask' )
parser.add_argument( '--source_mask_format', dest='source_mask_format', default=None, help='Source mask image format' )
parser.add_argument( '--target_image', dest='target_image', help='Target image' )
parser.add_argument( '--target_image_format', dest='target_image_format', help='Target image format' )
parser.add_argument( '--target_mask', dest='target_mask', default=None, help='Target mask' )
parser.add_argument( '--target_mask_format', dest='target_mask_format', default=None, help='Target mask image format' )
parser.add_argument( '--min_scale_def', dest='min_scale_def', type=int, help='Initial deformation' )
parser.add_argument( '--max_scale_def', dest='max_scale_def', type=int, help='Final deformation' )
parser.add_argument( '--max_subsamp_fact', dest='max_subsamp_fact', type=int, help='Image sub-sample factor' )
parser.add_argument( '--divergence_weight', dest='divergence_weight', type=float, help='Divergence weight' )
parser.add_argument( '--curl_weight', dest='curl_weight', type=float, help='Curl weight' )
parser.add_argument( '--image_weight', dest='image_weight', type=float, help='Image weight' )
parser.add_argument( '--consistency_weight', dest='consistency_weight', type=float, help='Consistency weight' )
parser.add_argument( '--landmarks_weight', dest='landmarks_weight', type=float, help='Landmarks weight' )
parser.add_argument( '--landmarks_file', dest='landmarks_file', default=None, help='Landmarks file' )
parser.add_argument( '--source_affine_file', dest='source_affine_file', default=None, help='Initial source affine matrix transformation' )
parser.add_argument( '--target_affine_file', dest='target_affine_file', default=None, help='Initial target affine matrix transformation' )
parser.add_argument( '--mono', dest='mono', default=False, help='Unidirectional registration (source to target)' )
parser.add_argument( '--source_trans_out', dest='source_trans_out', default=None, help='Direct source transformation matrix' )
parser.add_argument( '--target_trans_out', dest='target_trans_out', default=None, help='Inverse target transformation matrix' )
parser.add_argument( '--source_out', help='Output source image' )
parser.add_argument( '--source_out_datatype', help='Output registered source image format' )
parser.add_argument( '--target_out', default=None, help='Output target image' )
parser.add_argument( '--target_out_datatype', default=None, help='Output registered target image format' )
parser.add_argument( '--jython_script', dest='jython_script', help='Path to the Jython script' )

args = parser.parse_args()

if args.source_trans_out is not None and args.target_trans_out is not None:
    save_transformation = True
else:
    save_transformation = False

tmp_dir = imagej2_base_utils.get_temp_dir()
source_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.source_image, args.source_image_format )
tmp_source_out_tiff_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, 'tiff' )
tmp_source_out_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.source_out_datatype )
target_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.target_image, args.target_image_format )
if not args.mono:
    tmp_target_out_tiff_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, 'tiff' )
    tmp_target_out_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.target_out_datatype )
if args.source_mask is not None and args.target_mask is not None:
    tmp_source_mask_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.source_mask, args.source_mask_format )
    tmp_target_mask_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.target_mask, args.target_mask_format )
if save_transformation:
    # bUnwarpJ automatically names the transformation files based on the names
    # of the source and target image file names.  We've defined symlinks to 
    # temporary files with valid image extensions since ImageJ does not handle
    # the Galaxy "dataset.dat" file extensions.
    source_file_name = imagej2_base_utils.get_file_name_without_extension( tmp_source_out_tiff_path )
    tmp_source_out_transf_path = os.path.join( tmp_dir, '%s_transf.txt' % source_file_name )
    target_file_name = imagej2_base_utils.get_file_name_without_extension( tmp_target_out_tiff_path )
    tmp_target_out_transf_path = os.path.join( tmp_dir, '%s_transf.txt' % target_file_name )

# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open( tmp_out, 'wb' )
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open( tmp_err, 'wb' )

# Build the command line to align the two images.
cmd = imagej2_base_utils.get_base_cmd_bunwarpj( None )
if cmd is None:
    imagej2_base_utils.stop_err( "bUnwarpJ not found!" )
cmd += ' -align'
# Target is sent before source.
cmd += ' %s' % target_image_path
if args.target_mask is None:
    target_mask_str = ' NULL'
else:
    target_mask_str = ' %s' % tmp_target_mask_path
cmd += target_mask_str
cmd += ' %s' % source_image_path
if args.source_mask is None:
    source_mask_str = ' NULL'
else:
    source_mask_str = ' %s' % tmp_source_mask_path
cmd += source_mask_str
cmd += ' %d' % args.min_scale_def
cmd += ' %d' % args.max_scale_def
cmd += ' %d' % args.max_subsamp_fact
cmd += ' %.1f' % args.divergence_weight
cmd += ' %.1f' % args.curl_weight
cmd += ' %.1f' % args.image_weight
cmd += ' %.1f' % args.consistency_weight
# Source is produced before target.
cmd += ' %s' % tmp_source_out_tiff_path
if not args.mono:
    cmd += ' %s' % tmp_target_out_tiff_path
if args.landmarks_file is not None:
    # We have to create a temporary file with a .txt extension here so that
    # bUnwarpJ will not ignore the Galaxy "dataset.dat" file.
    tmp_landmarks_file_path = imagej2_base_utils.get_input_image_path( tmp_dir,
                                                                       args.landmarks_file,
                                                                       'txt' )
    cmd += ' -landmarks'
    cmd += ' %.1f' % args.landmarks_weight
    cmd += ' %s' % tmp_landmarks_file_path
if args.source_affine_file is not None and args.target_affine_file is not None:
    # Target is sent before source.
    cmd += ' -affine'
    cmd += ' %s' % args.target_affine_file
    cmd += ' %s' % args.source_affine_file
if args.mono:
    cmd += ' -mono'
if save_transformation:
    cmd += ' -save_transformation'

# Align the two images using bUnwarpJ.
proc = subprocess.Popen( args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True )
rc = proc.wait()
if rc != 0:
    error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
    imagej2_base_utils.stop_err( error_message )

# bUnwarpJ produces tiff image stacks consisting of 3 slices which can be viewed in ImageJ.
# The 3 slices are:: 1) the registered image, 2) the target image and 3) the black/white
# warp image.  Galaxy supports only single-layered images, so we'll convert the images so they
# can be viewed in Galaxy.

# Define command response buffers.
tmp_out = tempfile.NamedTemporaryFile().name
tmp_stdout = open( tmp_out, 'wb' )
tmp_err = tempfile.NamedTemporaryFile().name
tmp_stderr = open( tmp_err, 'wb' )

# Build the command line to handle the multi-slice tiff images.
cmd = imagej2_base_utils.get_base_command_imagej2( None, jython_script=args.jython_script )
if cmd is None:
    imagej2_base_utils.stop_err( "ImageJ not found!" )
if args.mono:
    # bUnwarpJ will produce only a registered source image.
    cmd += ' %s %s %s %s' % ( tmp_source_out_tiff_path,
                              args.source_out_datatype,
                              tmp_source_out_path,
                              args.mono )
else:
    # bUnwarpJ will produce registered source and target images.
    cmd += ' %s %s %s %s %s %s %s' % ( tmp_source_out_tiff_path,
                                       args.source_out_datatype,
                                       tmp_source_out_path,
                                       tmp_target_out_tiff_path,
                                       args.target_out_datatype,
                                       tmp_target_out_path,
                                       args.mono )

# Merge the multi-slice tiff layers into an image that can be viewed in Galaxy.
proc = subprocess.Popen( args=cmd, stderr=tmp_stderr, stdout=tmp_stdout, shell=True )
rc = proc.wait()
if rc != 0:
    error_message = imagej2_base_utils.get_stderr_exception( tmp_err, tmp_stderr, tmp_out, tmp_stdout )
    imagej2_base_utils.stop_err( error_message )

# Save the Registered Source Image to the output dataset.
shutil.move( tmp_source_out_path, args.source_out )
if not args.mono:
    # Move the Registered Target Image to the output dataset.
    shutil.move( tmp_target_out_path, args.target_out )

# If requested, save matrix transformations as additional datasets.
if save_transformation:
    shutil.move( tmp_source_out_transf_path, args.source_trans_out )
    if not args.mono:
        shutil.move( tmp_target_out_transf_path, args.target_trans_out )

imagej2_base_utils.cleanup_before_exit( tmp_dir )
