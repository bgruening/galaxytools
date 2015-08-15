import jython_utils
import sys
from ij import IJ
from ij import ImagePlus
from ij.plugin.filter import Analyzer
from ij.plugin.filter import EDM

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -8 ]
input = sys.argv[ -7 ]
iterations = int( sys.argv[ -6 ] )
count = int( sys.argv[ -5 ] )
black_background = jython_utils.asbool( sys.argv[ -4 ] )
pad_edges_when_eroding = jython_utils.asbool( sys.argv[ -3 ] )
tmp_output_path = sys.argv[ -2 ]
output_datatype = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()

try:
    # Set binary options.
    options = jython_utils.get_binary_options( black_background=black_background,
                                               iterations=iterations,
                                               count=count,
                                               pad_edges_when_eroding=pad_edges_when_eroding )
    IJ.run( input_image_plus_copy, "Options...", options )

    # Convert image to binary if necessary.
    if not image_processor_copy.isBinary():
        # Convert the image to binary grayscale.
        IJ.run( input_image_plus_copy, "Make Binary", "" )

    # Run the command.
    IJ.run( input_image_plus_copy, "Distance Map", "" )
    # Save the ImagePlus object as a new image.
    IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
