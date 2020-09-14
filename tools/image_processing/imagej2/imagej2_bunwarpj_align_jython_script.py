import sys

from ij import IJ


# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
mono = sys.argv[-1] == "true"

if mono:
    # bUnwarpJ has been called with the -mono param.
    source_tiff_path = sys.argv[-4]
    source_datatype = sys.argv[-3]
    source_path = sys.argv[-2]
else:
    source_tiff_path = sys.argv[-7]
    source_datatype = sys.argv[-6]
    source_path = sys.argv[-5]
    target_tiff_path = sys.argv[-4]
    target_datatype = sys.argv[-3]
    target_path = sys.argv[-2]

# Save the Registered Source Image.
registered_source_image = IJ.openImage(source_tiff_path)
IJ.saveAs(registered_source_image, source_datatype, source_path)

if not mono:
    # Save the Registered Target Image.
    registered_target_image = IJ.openImage(target_tiff_path)
    IJ.saveAs(registered_target_image, target_datatype, target_path)
