import sys

from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
input_file = sys.argv[-8]
image_datatype = sys.argv[-7]
filter = sys.argv[-6]
radius = sys.argv[-5]
mask = sys.argv[-4]
light_background = sys.argv[-3]
dont_substract = sys.argv[-2]
tmp_output_path = sys.argv[-1]

# Open the input image file.
image_plus = IJ.openImage(input_file)
# Create an ImagePlus object for the image.
image_plus_copy = image_plus.duplicate()

# Perform the analysis on the ImagePlus object.
try:
    if filter == "gaussian_blur":
        IJ.run(image_plus_copy, "Gaussian Blur...", "sigma=%s" % radius)
    elif filter in ["median", "mean", "minimum", "maximum", "variance"]:
        IJ.run(image_plus_copy, "%s..." % filter.title(), "radius=%s" % radius)
    elif filter == "unsharp_mask":
        IJ.run(image_plus_copy, "Unsharp Mask...", "radius=%s mask=%s" % (radius, mask))
    elif filter == "top_hat":
        print("radius=%s %s %s" % (radius, light_background, dont_substract.replace('dont', "don't")))
        IJ.run(image_plus_copy, "Top Hat...", "radius=%s %s %s" % (radius, light_background, dont_substract.replace('dont', "don't")))
except Exception as e:
    # This is due to some operations like gaussian_blur which block the script
    print(e)
    exit(1)
# Save the ImagePlus object as a new image.
IJ.saveAs(image_plus_copy, image_datatype, tmp_output_path)
# This is due to some operations like gaussian_blur which block the script
exit(0)
