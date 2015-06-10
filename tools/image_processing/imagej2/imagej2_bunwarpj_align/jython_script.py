import sys
import imagej2_base_utils
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

def convert_before_saving_as_tiff( registered_image ):
    # The bUnwarpJ macro produces tiff image stacks consisting of 3
    # slices which can be viewed in ImageJ.  The 3 slices are::
    # 1) the registered image, 2) the target image and 3) the black/white
    # warp image.  Galaxy supports only single-layered images, so we have
    # to convert the image to something other than tiff so that slices are
    # eliminated and we can then convert back to tiff for saving.  This is
    # kind of a hack - there is probably a way to do this without converting
    # twice.  I messed around with the the ImageJ getSlice() method, but
    # couldn't get it to produce what is needed here.
    tmp_dir = imagej2_base_utils.get_temp_dir()
    tmp_out_png_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, 'png' )
    IJ.saveAs( registered_image, 'png', tmp_out_png_path )
    return IJ.openImage( tmp_out_png_path )

# Save the Registered Source Image.
registered_source_image = IJ.openImage( source_tiff_path )
if source_datatype == 'tiff':
    registered_source_image = convert_before_saving_as_tiff( registered_source_image )
IJ.saveAs( registered_source_image, source_datatype, source_path )

if not mono:
    # Save the Registered Target Image.
    registered_target_image = IJ.openImage( target_tiff_path )
    if target_datatype == 'tiff':
        registered_target_image = convert_before_saving_as_tiff( registered_target_image )
    IJ.saveAs( registered_target_image, target_datatype, target_path )
