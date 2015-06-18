import sys
import jython_utils
from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.

source_tiff_path = sys.argv[ -3 ]
source_datatype = sys.argv[ -2 ]
source_path = sys.argv[ -1 ]

# Save the Registered Source Image.
registered_source_image = IJ.openImage( source_tiff_path )
if source_datatype == 'tiff':
    registered_source_image = jython_utils.convert_before_saving_as_tiff( registered_source_image )
IJ.saveAs( registered_source_image, source_datatype, source_path )
