import sys

from ij import IJ, ImagePlus

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
input_file = sys.argv[-8]
image_datatype = sys.argv[-7]
noise = sys.argv[-6]
standard_deviation = sys.argv[-5]
radius = sys.argv[-4]
threshold = sys.argv[-3]
which_outliers = sys.argv[-2]
tmp_output_path = sys.argv[-1]

# Open the input image file.
image_plus = IJ.openImage(input_file)
image_type = image_plus.getType()
is32BITS_GREY = image_type == ImagePlus.GRAY32
# Create an ImagePlus object for the image.
image_plus_copy = image_plus.duplicate()

# Perform the analysis on the ImagePlus object.
try:
    if noise == "add_noise":
        IJ.run(image_plus_copy, "Add Noise", "")
    elif noise == "add_specified_noise":
        IJ.run(image_plus_copy, "Add Specified Noise...", "standard=%s" % standard_deviation)
    elif noise == "salt_and_pepper":
        IJ.run(image_plus_copy, "Salt and Pepper", "")
    elif noise == "despeckle":
        IJ.run(image_plus_copy, "Despeckle", "")
    elif noise == "remove_outliers":
        IJ.run(
            image_plus_copy,
            "Remove Outliers...",
            "radius=%s threshold=%s which=%s" % (radius, threshold, which_outliers)
        )
    elif noise == "remove_nans":
        if is32BITS_GREY:
            IJ.run(image_plus_copy, "Remove NaNs...", "")
        else:
            raise Exception("Remove NaNs can only be applied to 32bits grey images.")
    elif noise == "rof_denoise":
        if is32BITS_GREY:
            IJ.run(image_plus_copy, "ROF Denoise", "")
        else:
            raise Exception("ROF Denoise can only be applied to 32bits grey images.")
except Exception as e:
    # This is due to some operations like remove_outliers and despeckle which block the script
    print(e)
    exit(1)
# Save the ImagePlus object as a new image.
IJ.saveAs(image_plus_copy, image_datatype, tmp_output_path)
# This is due to some operations like remove_outliers and despeckle which block the script
exit(0)
