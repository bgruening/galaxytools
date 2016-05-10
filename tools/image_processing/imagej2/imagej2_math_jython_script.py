import jython_utils
import sys
from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -8 ]
input = sys.argv[ -7 ]
operation = sys.argv[ -6 ]
expression = sys.argv[ -5 ]
if sys.argv[ -4 ] in [ None, 'None' ]:
    bin_constant = None
else:
    bin_constant = int( sys.argv[ -4 ] )
if sys.argv[ -3 ] in [ None, 'None' ]:
    float_constant = None
else:
    float_constant = float( sys.argv[ -3 ] )
tmp_output_path = sys.argv[ -2 ]
output_datatype = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()
bit_depth = image_processor_copy.getBitDepth()

try:
    if operation.find( '_' ) > 0:
        # Square_Root.
        new_operation = operation.replace( '_', ' ' )
    elif operation in [ 'Square', 'Log', 'Exp', 'Abs', 'Reciprocal' ]:
        # Unfortunately some ImageJ commands require a "..." ending
        # while others do not.  There seems to be no pattern.
        new_operation = '%s' % operation
    else:
        new_operation = '%s...' % operation

    if operation == 'Macro':
        # Apply the macro code to the image via a call to it's
        # ImageProcessor since this option does not work using
        # the IJ.run() method.
        new_expression = expression.lstrip( '"' ).rstrip( '"' )
        options = 'code=%s' % new_expression
        image_processor_copy.applyMacro( new_expression )
    elif operation == 'Min':
        # Min does not work without using the ImageProcessor.
        image_processor_copy.min( float_constant )
    elif operation == 'Max':
        # Max does not work without using the ImageProcessor.
        image_processor_copy.max( float_constant )
    elif operation == 'Abs':
        if bit_depth not in [ 16, 32 ]:
            # Convert the image to 32-bit.
            IJ.run( input_image_plus_copy, "32-bit", "" )
            IJ.run( input_image_plus_copy, new_operation, "" )
    elif operation == 'Reciprocal':
        if bit_depth != 32:
            # Convert the image to 32 bit.
            IJ.run( input_image_plus_copy, "32-bit", "" )
            IJ.run( input_image_plus_copy, new_operation, "" )
    else:
        if operation in [ 'AND', 'OR', 'XOR' ]:
            # Value is a binary number.
            options = 'value=%d' % bin_constant
        elif operation in [ 'Log', 'Exp', 'Square', 'Square_Root' ]:
            # No constant value.
            options = ''
        else:
            # Value is a floating point number.
            options = 'value=%.3f' % float_constant
        IJ.run( input_image_plus_copy, "%s" % new_operation, "%s" % options )
    # Save the ImagePlus object as a new image.
    IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
