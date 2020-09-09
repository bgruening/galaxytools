import os
import sys
import tempfile

from ij import IJ


# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.

source_tiff_path = sys.argv[-3]
source_datatype = sys.argv[-2]
source_path = sys.argv[-1]

# Save the Registered Source Image.
registered_source_image = IJ.openImage(source_tiff_path)
if source_datatype == 'tiff':
    # The bUnwarpJ plug-in produces TIFF image stacks consisting of 3
    # slices which can be viewed in ImageJ.  The 3 slices are: 1) the
    # registered image, 2) the target image and 3) the black/white warp
    # image.  When running bUnwarpJ from the command line (as these
    # Galaxy wrappers do) the initial call to IJ.openImage() (to open the
    # registered source and target images produced by bUnwarpJ) in the
    # tool's jython_script.py returns an ImagePlus object with a single
    # slice which is the "generally undesired" slice 3 discussed above.
    # However, a call to IJ.saveAs() will convert the single-slice TIFF
    # into a 3-slice TIFF image stack (as described above) if the selected
    # format for saving is TIFF.  Galaxy supports only single-layered
    # images, so to work around this behavior, we have to convert the
    # image to something other than TIFF so that slices are eliminated.
    # We can then convert back to TIFF for saving.  There might be a way
    # to do this without converting twice, but no idea.
    tmp_dir = tempfile.mkdtemp(prefix='tmp-imagej-')
    fd, tmp_out_png_path = tempfile.mkstemp(suffix='png', dir=tmp_dir)
    os.close(fd)
    IJ.saveAs(registered_source_image, 'png', tmp_out_png_path)
    registered_source_image = IJ.openImage(tmp_out_png_path)
IJ.saveAs(registered_source_image, source_datatype, source_path)
