import sys
import jython_utils
from ij import ImagePlus, IJ
from ij.plugin.filter import Analyzer, MaximumFinder
from ij.process import ImageProcessor
from jarray import array

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -10 ]
input = sys.argv[ -9 ]
scale_when_converting = jython_utils.asbool( sys.argv[ -8 ] )
weighted_rgb_conversions = jython_utils.asbool( sys.argv[ -7 ] )
noise_tolerance = int( sys.argv[ -6 ] )
output_type = sys.argv[ -5 ]
exclude_edge_maxima = jython_utils.asbool( sys.argv[ -4 ] )
light_background = jython_utils.asbool( sys.argv[ -3 ] )
tmp_output_path = sys.argv[ -2 ]
output_datatype = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()
bit_depth = image_processor_copy.getBitDepth()
analyzer = Analyzer( input_image_plus_copy )

try:
    # Set the conversion options.
    options = []
    # The following 2 options are applicable only to RGB images.
    if bit_depth == 24:
        if scale_when_converting:
            option.append( "scale" )
        if weighted_rgb_conversions:
            options.append( "weighted" )
    # Perform conversion - must happen even if no options are set.
    IJ.run( input_image_plus_copy, "Conversions...", "%s" % " ".join( options ) )
    if output_type in [ 'List', 'Count' ]:
        # W're  generating a tabular file for the output.
        # Set the Find Maxima options.
        options = [ 'noise=%d' % noise_tolerance ]
        if output_type.find( '_' ) > 0:
            output_type_str = 'output=[%s]' % output_type.replace( '_', ' ' )
        else:
            output_type_str = 'output=%s' % output_type
        options.append( output_type_str )
        if exclude_edge_maxima:
            options.append( 'exclude' )
        if light_background:
            options.append( 'light' )
        # Run the command.
        IJ.run( input_image_plus_copy, "Find Maxima...", "%s" % " ".join( options ) )
        results_table = analyzer.getResultsTable()
        results_table.saveAs( tmp_output_path )
    else:
        # Find the maxima of an image (does not find minima).
        # LIMITATIONS: With output_type=Segmented_Particles
        # (watershed segmentation), some segmentation lines
        # may be improperly placed if local maxima are suppressed
        # by the tolerance.
        mf = MaximumFinder()
        if output_type == 'Single_Points':
            output_type_param = mf.SINGLE_POINTS
        elif output_type == 'Maxima_Within_Tolerance':
            output_type_param = mf.IN_TOLERANCE
        elif output_type == 'Segmented_Particles':
            output_type_param = mf.SEGMENTED
        elif output_type == 'List':
            output_type_param = mf.LIST
        elif output_type == 'Count':
            output_type_param = mf.COUNT
        # Get a new byteProcessor with a normal (uninverted) LUT where
        # the marked points are set to 255 (Background 0). Pixels outside
        # of the roi of the input image_processor_copy are not set. No
        # output image is created for output types POINT_SELECTION, LIST
        # and COUNT.  In these cases findMaxima returns null.
        byte_processor = mf.findMaxima( image_processor_copy,
                                        noise_tolerance,
                                        ImageProcessor.NO_THRESHOLD,
                                        output_type_param,
                                        exclude_edge_maxima,
                                        False )
        # Invert the image or ROI.
        byte_processor.invert()
        if output_type == 'Segmented_Particles' and not light_background:
            # Invert the values in this image's LUT (indexed color model).
            byte_processor.invertLut()
        image_plus = ImagePlus( "output", byte_processor )
        IJ.saveAs( image_plus, output_datatype, tmp_output_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
