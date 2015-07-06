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
if image_type in VALID_IMAGE_TYPES:
    # Create a copy of the image.
    input_image_plus_copy = input_image_plus.createImagePlus()
    image_processor_copy = input_image_plus.getProcessor().duplicate()
    input_image_plus_copy.setProcessor( "iCopy", image_processor_copy )
    # Run the command.
    IJ.run( input_image_plus_copy, "Skeletonize (2D/3D)", "" )
    # Save the ImagePlus object as a new image.
    IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
else:
    # When Galaxy metadata for images is enhanced to include information like this,
    # we'll be able to write tool validators rather than having to stop the job in
    # an error state.
    display_image_type = jython_utils.get_display_image_type( image_type )
    if display_image_type is None:
        display_image_type = image_type
    msg = "Skeletonize3D requires an 8-bit grayscale image, the selected image is %s." % display_image_type
    jython_utils.handle_error( error_log, msg )
