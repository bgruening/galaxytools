#!/usr/bin/env python
# Galaxy wrapper for use with Imagej2 via bioformats and javabridge by Greg Von Kuster
"""
A wrapper script for running ImageJ2 commands via bioformats and javabridge.
"""
import optparse
import imagej2_utils

if __name__=="__main__":
    # Parse Command Line.
    parser = optparse.OptionParser()
    parser.add_option( '-i', '--in_fname', dest='in_fname', action='store', type='string', default=None )
    parser.add_option( '-d', '--input_datatype', dest='input_datatype', action='store', type="string", default='bmp' )
    parser.add_option( '-t', '--max_heap_size_type', dest='max_heap_size_type', action='store', type="string", default='default' )
    parser.add_option( '-m', '--max_heap_size', dest='max_heap_size', action='store', type="int", default=0, help='Maximum size of the memory allocation pool used by the JVM.' )
    parser.add_option( '-q', '--output_datatype', dest='output_datatype', action='store', type="string", default='jpg' )
    parser.add_option( '-o', '--out_fname', dest='out_fname', action='store', type='string', default=None )
    ( options, args ) = parser.parse_args()
    # Set the size of the memory allocation pool used by the JVM.
    max_heap_size = imagej2_utils.get_max_heap_size_value( options.max_heap_size_type, options.max_heap_size )
    # Start the JVM via the Javabridge.
    imagej2_utils.start_vm( args=None, class_path=None, max_heap_size=max_heap_size, run_headless=True )
    try:
        tmp_dir = imagej2_utils.get_temp_dir()
        in_image_path = imagej2_utils.get_input_image_path( tmp_dir, options.in_fname, options.input_datatype )
        # Load the input image.
        image, scale = imagej2_utils.load_image( in_image_path, rescale=False, wants_max_intensity=True )
        # Write the output image.
        out_image_path = imagej2_utils.get_temporary_image_path( tmp_dir, options.output_datatype )
        imagej2_utils.write_image( image_path=out_image_path, pixels=image, pixel_type=str( image.dtype ), move_to=options.out_fname )
    except Exception, e:
        imagej2_utils.stop_err( str( e ) )
    finally:
        imagej2_utils.kill_vm()
        imagej2_utils.cleanup_before_exit( tmp_dir )
