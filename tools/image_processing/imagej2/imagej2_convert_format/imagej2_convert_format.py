#!/usr/bin/env python
import argparse
import imagej2_base_utils
import imagej2_utils

if __name__=="__main__":
    # Parse Command Line.
    parser = argparse.ArgumentParser()
    parser.add_argument( '--in_fname', dest='in_fname', help='Path to the input file' )
    parser.add_argument( '--input_datatype', dest='input_datatype', help='Input image datatype' )
    parser.add_argument( '--max_heap_size_type', dest='max_heap_size_type', help='Type (default or megabytes) of max_heap_size value' )
    parser.add_argument( '--max_heap_size', dest='max_heap_size', help='Maximum size of the memory allocation pool used by the JVM.' )
    parser.add_argument( '--output_datatype', dest='output_datatype', help='Output image datatype' )
    parser.add_argument( '--out_fname', help='Path to the output file' )
    args = parser.parse_args()
    # Set the size of the memory allocation pool used by the JVM.
    max_heap_size = imagej2_base_utils.get_max_heap_size_value( args.max_heap_size_type, args.max_heap_size )
    # Start the JVM via the Javabridge.
    imagej2_utils.start_vm( args=None, class_path=None, max_heap_size=max_heap_size, run_headless=True )
    try:
        tmp_dir = imagej2_base_utils.get_temp_dir()
        in_image_path = imagej2_base_utils.get_input_image_path( tmp_dir, args.in_fname, args.input_datatype )
        # Load the input image.
        image, scale = imagej2_utils.load_image( in_image_path, rescale=False, wants_max_intensity=True )
        # Write the output image.
        out_image_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, args.output_datatype )
        imagej2_utils.write_image( image_path=out_image_path, pixels=image, pixel_type=str( image.dtype ), move_to=args.out_fname )
    except Exception, e:
        imagej2_base_utils.stop_err( str( e ) )
    finally:
        imagej2_utils.kill_vm()
        imagej2_base_utils.cleanup_before_exit( tmp_dir )
