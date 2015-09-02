import jython_utils
import sys
from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -7 ]
input = sys.argv[ -6 ]
equalize_histogram = jython_utils.asbool( sys.argv[ -5 ] )
saturated_pixels = sys.argv[ -4 ]
normalize = jython_utils.asbool( sys.argv[ -3 ] )
tmp_output_path = sys.argv[ -2 ]
output_datatype = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()
bit_depth = image_processor_copy.getBitDepth()

# Set the options
options = []
# If equalize_histogram, saturated_pixels and normalize are ignored.
if equalize_histogram:
    options.append( 'equalize' )
else:
    if saturated_pixels not in [ None, 'None' ]:
        # Fiji allows only a single decimal place for this value.
        options.append( 'saturated=%.3f' % float( saturated_pixels ) )
    # Normalization of RGB images is not supported.
    if bit_depth != 24 and normalize:
        options.append( 'normalize' )
try:
    # Run the command.
    options = "%s" % ' '.join( options )
    IJ.run( input_image_plus_copy, "Enhance Contrast...", options )
    # Save the ImagePlus object as a new image.
    IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
