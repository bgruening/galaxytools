import sys
import jython_utils
from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.

if sys.argv[ -1 ].lower() in [ 'true' ]:
    mono = True
else:
    mono = False

if mono:
    # bUnwarpJ has been called with the -mono param.
    source_tiff_path = sys.argv[ -4 ]
    source_datatype = sys.argv[ -3 ]
    source_path = sys.argv[ -2 ]
else:
    source_tiff_path = sys.argv[ -7 ]
    source_datatype = sys.argv[ -6 ]
    source_path = sys.argv[ -5 ]
    target_tiff_path = sys.argv[ -4 ]
    target_datatype = sys.argv[ -3 ]
    target_path = sys.argv[ -2 ]

# Save the Registered Source Image.
registered_source_image = IJ.openImage( source_tiff_path )
if source_datatype == 'tiff':
    registered_source_image = jython_utils.convert_before_saving_as_tiff( registered_source_image )
IJ.saveAs( registered_source_image, source_datatype, source_path )

if not mono:
    # Save the Registered Target Image.
    registered_target_image = IJ.openImage( target_tiff_path )
    if target_datatype == 'tiff':
        registered_target_image = jython_utils.convert_before_saving_as_tiff( registered_target_image )
    IJ.saveAs( registered_target_image, target_datatype, target_path )
