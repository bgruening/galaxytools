import imagej2_base_utils
from ij import IJ

def convert_before_saving_as_tiff( image_plus ):
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
    # to do this without converting twice, but I spent a lot of time looking
    # and I have yet to discover it.
    tmp_dir = imagej2_base_utils.get_temp_dir()
    tmp_out_png_path = imagej2_base_utils.get_temporary_image_path( tmp_dir, 'png' )
    IJ.saveAs( image_plus, 'png', tmp_out_png_path )
    return IJ.openImage( tmp_out_png_path )
