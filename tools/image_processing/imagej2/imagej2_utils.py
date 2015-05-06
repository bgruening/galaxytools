import bioformats
import javabridge
import os
import shutil
import subprocess
import sys
import tempfile

def cleanup_before_exit( tmp_dir ):
    """
    Remove temporary files and directories prior to tool exit.
    """
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )

def get_file_extension( image_format ):
    """
    Return a valid bioformats file extension based on the received
    value of image_format( e.g., "gif" is returned as ".gif".
    """
    return '.%s' % image_format

def get_temporary_image_path( tmp_dir, image_format ):
    """
    Return the path to a temporary file with a valid image format
    file extension that can be used with bioformats.
    """
    file_extension = get_file_extension( image_format )
    return get_tempfilename( tmp_dir, file_extension )

def get_input_image_path( tmp_dir, input_file, image_format ):
    """
    Bioformats uses file extensions (e.g., .job, .gif, etc)
    when reading and writing image files, so the Galaxy dataset
    naming convention of setting all file extensions as .dat
    must be handled.
    """
    image_path = get_temporary_image_path( tmp_dir, image_format )
    # Remove the file so we can create a symlink.
    os.remove( image_path )
    os.symlink( input_file, image_path )
    return image_path

def get_java_class_path():
    """
    Return the Java class path for use when starting the JVM via Javabridge.
    This function expects the following environment variable settings:
    FIJI_JAR_DIR, FIJI_PLUGIN_DIR, FIJI_OSX_JAVA3D_DIR.
    """
    class_path_list = []
    # Handle javabridge.JARS setting.
    for jar_file in javabridge.JARS:
        class_path_list.append( jar_file )
    fiji_jar_dir = os.environ.get( 'FIJI_JAR_DIR', None )
    if fiji_jar_dir is not None and os.path.isdir( fiji_jar_dir ):
        for filename in os.listdir( fiji_jar_dir ):
            if filename.endswith( '.jar' ):
                class_path_list.append( os.path.join( fiji_jar_dir, filename ) )
        fiji_bioformats_jar_dir = os.path.join( fiji_jar_dir, 'bio-formats' )
        if os.path.isdir( fiji_bioformats_jar_dir ):
            for filename in os.listdir( fiji_bioformats_jar_dir ):
                if filename.endswith( '.jar' ):
                    class_path_list.append( os.path.join( fiji_bioformats_jar_dir, filename ) )
    fiji_plugin_dir = os.environ.get( 'FIJI_PLUGIN_DIR', None )
    if fiji_plugin_dir is not None and os.path.isdir( fiji_plugin_dir ):
        for filename in os.listdir( fiji_plugin_dir ):
            if filename.endswith( '.jar' ):
                class_path_list.append( os.path.join( fiji_plugin_dir, filename ) )
    fiji_osx_java3d_dir = os.environ.get( 'FIJI_OSX_JAVA3D_DIR', None )
    if fiji_osx_java3d_dir is not None and os.path.isdir( fiji_osx_java3d_dir ):
        for filename in os.listdir( fiji_osx_java3d_dir ):
            if filename.endswith( '.jar' ):
                class_path_list.append( os.path.join( fiji_osx_java3d_dir, filename ) )
    return class_path_list

def get_max_heap_size_value( max_heap_size_type, max_heap_size ):
    """
    Return a string that can be used by the javabridge to set the size
    of the memory allocation pool used by the JVM.  The value must be
    determined to be a multiple of 1024 or it will be ignored.
    """
    if max_heap_size_type == 'default':
        return None
    if max_heap_size_type == 'megabytes':
        if max_heap_size % 1024 not in [ 0, 256, 512 ]:
            return None
        return '%sm' % str( max_heap_size )

def get_temp_dir( prefix='tmp-imagej-' ):
    """
    Return a temporary directory.
    """
    return tempfile.mkdtemp( prefix=prefix )

def get_tempfilename( dir, suffix ):
    """
    Return a temporary file name.
    """
    fd, name = tempfile.mkstemp( suffix, dir )
    os.close( fd )
    return name

def kill_vm():
    javabridge.detach()
    javabridge.kill_vm()

def load_image( image_path, c=None, z=0, t=0, series=None, index=None,
                rescale=True, wants_max_intensity=False, channel_names=None ):
    """
    Load in image from a file using bioformats, always returning a tuple of
    (image, scale) where the image is a numpy.ndarray object, and the scale
    is an integer or None.
    """
    result = bioformats.load_image( path=image_path,
                                    c=c,
                                    z=z,
                                    t=t,
                                    series=series,
                                    index=index,
                                    rescale=rescale,
                                    wants_max_intensity=wants_max_intensity,
                                    channel_names=channel_names )
    if wants_max_intensity:
        # The value of result is a tuple: ( image, scale )
        image, scale = result
    else:
        image = result
        scale = None
     # The image is a numpy.ndarray object, and the scale is an integer or None.
    return image, scale

def start_vm( args=None, class_path=None, max_heap_size=None, run_headless=False ):
    """
    Start the JVM via Javabridge.
    """
    if class_path is None:
        class_path = get_java_class_path()
    # Start the JVM.
    javabridge.start_vm( args=args, class_path=class_path, max_heap_size=max_heap_size, run_headless=run_headless )
    # Suppress Java logging to stdout.
    java_stack = javabridge.make_instance( 'java/io/ByteArrayOutputStream', "()V" )
    java_stack_ps = javabridge.make_instance( 'java/io/PrintStream', "(Ljava/io/OutputStream;)V", java_stack )
    javabridge.static_call( 'Ljava/lang/System;', "setErr", '(Ljava/io/PrintStream;)V', java_stack_ps )
    java_out = javabridge.make_instance( 'java/io/ByteArrayOutputStream', "()V" )
    java_out_ps = javabridge.make_instance( 'java/io/PrintStream', "(Ljava/io/OutputStream;)V", java_out )
    javabridge.static_call( 'Ljava/lang/System;', "setOut", '(Ljava/io/PrintStream;)V', java_out_ps )
    javabridge.attach()

def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit( 1 )

def write_image( image_path, pixels, pixel_type, c=0, z=0, t=0,
                 size_c=1, size_z=1, size_t=1, channel_names=None,
                 move_to=None ):
    """
    Write an image file using bioformats.  Bioformats uses file extensions
    (e.g., .job, .gif, etc) when reading and writing image files, so the
    Galaxy dataset naming convention of setting all file extensions as .dat
    is handled by setting the move_to parameter to the value of the Galaxy
    dataset path.
    """
    bioformats.write_image( pathname=image_path,
                            pixels=pixels,
                            pixel_type=pixel_type,
                            c=c,
                            z=z,
                            t=t,
                            size_c=size_c,
                            size_z=size_z,
                            size_t=size_t,
                            channel_names=channel_names )
    if move_to is not None and os.path.exists( move_to ):
        shutil.move( image_path, os.path.abspath( move_to ) )
