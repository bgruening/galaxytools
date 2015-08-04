import jython_utils
import sys
from ij import IJ
from ij import ImagePlus

VALID_IMAGE_TYPES = [ ImagePlus.GRAY8 ]

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -4 ]
input = sys.argv[ -3 ]
tmp_output_path = sys.argv[ -2 ]
output_datatype = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )
bit_depth = input_image_plus.getBitDepth()
image_type = input_image_plus.getType()

# Create a copy of the image.
input_image_plus_copy = input_image_plus.createImagePlus()
image_processor_copy = input_image_plus.getProcessor().duplicate()
input_image_plus_copy.setProcessor( "iCopy", image_processor_copy )

try:
    if image_type not in VALID_IMAGE_TYPES:
        # Convert the image to binary grayscale.
        IJ.run( input_image_plus_copy, "Make Binary","iterations=1 count=1 edm=Overwrite do=Nothing" )
    # Run the command.
    IJ.run( input_image_plus_copy, "Skeletonize (2D/3D)", "" )
    # Save the ImagePlus object as a new image.
    IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
