import bioformats
import javabridge
import os
import shutil

FIJI_JAR_DIR = os.environ.get( 'FIJI_JAR_DIR', None )
FIJI_OSX_JAVA3D_DIR = os.environ.get( 'FIJI_OSX_JAVA3D_DIR', None )
FIJI_PLUGIN_DIR = os.environ.get( 'FIJI_PLUGIN_DIR', None )
FIJI_ROOT_DIR = os.environ.get( 'FIJI_ROOT_DIR', None )

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
    if FIJI_JAR_DIR is not None and os.path.isdir( FIJI_JAR_DIR ):
        for filename in os.listdir( FIJI_JAR_DIR ):
            if filename.endswith( '.jar' ):
                class_path_list.append( os.path.join( FIJI_JAR_DIR, filename ) )
        fiji_bioformats_jar_dir = os.path.join( FIJI_JAR_DIR, 'bio-formats' )
        if os.path.isdir( fiji_bioformats_jar_dir ):
            for filename in os.listdir( fiji_bioformats_jar_dir ):
                if filename.endswith( '.jar' ):
                    class_path_list.append( os.path.join( fiji_bioformats_jar_dir, filename ) )
    if FIJI_PLUGIN_DIR is not None and os.path.isdir( FIJI_PLUGIN_DIR ):
        for filename in os.listdir( FIJI_PLUGIN_DIR ):
            if filename.endswith( '.jar' ):
                class_path_list.append( os.path.join( FIJI_PLUGIN_DIR, filename ) )
    if FIJI_OSX_JAVA3D_DIR is not None and os.path.isdir( FIJI_OSX_JAVA3D_DIR ):
        for filename in os.listdir( FIJI_OSX_JAVA3D_DIR ):
            if filename.endswith( '.jar' ):
                class_path_list.append( os.path.join( FIJI_OSX_JAVA3D_DIR, filename ) )
    return class_path_list

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
