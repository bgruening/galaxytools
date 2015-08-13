import jython_utils
import sys
from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -8 ]
input = sys.argv[ -7 ]
iterations = sys.argv[ -6 ]
count = sys.argv[ -5 ]
black_background = sys.argv[ -4 ]
pad_edges_when_eroding = sys.argv[ -3 ]
tmp_output_path = sys.argv[ -2 ]
output_datatype = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()

try:
    options = [ "iterations=&iterations" "count=&count" ]
    if black_background == "yes":
        options.append( "black" )
    if pad_edges_when_eroding == "yes":
        options.append( "pad" )
    # The following have no affect within this tool, but are necessary.
    options.append( "edm=Overwrite" )
    options.append( "do=Nothing" )
    # Run the command.
    IJ.run( input_image_plus_copy, "Make Binary", " ".join( options ) )
    # Save the ImagePlus object as a new image.
    IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
