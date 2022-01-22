import sys

from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[-10]
input_file = sys.argv[-9]
threshold_min = float(sys.argv[-8])
threshold_max = float(sys.argv[-7])
method = sys.argv[-6]
display = sys.argv[-5]
black_background = sys.argv[-4] == "yes"
# TODO: this is not being used.
stack_histogram = sys.argv[-3] == "yes"
output_filename = sys.argv[-2]
output_datatype = sys.argv[-1]

# Open the input image file.
input_image_plus = IJ.openImage(input_file)

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()

# Convert image to binary if necessary.
if not image_processor_copy.isBinary():
    # Convert the image to binary grayscale.
    IJ.run(
        input_image_plus_copy,
        "Make Binary",
        "iterations=1 count=1 edm=Overwrite do=Nothing",
    )
# Set the options.
if black_background:
    method_str = "%s dark" % method
else:
    method_str = method
IJ.setAutoThreshold(input_image_plus_copy, method_str)
if display == "red":
    display_mode = "Red"
elif display == "bw":
    display_mode = "Black & White"
elif display == "over_under":
    display_mode = "Over/Under"
IJ.setThreshold(input_image_plus_copy, threshold_min, threshold_max, display_mode)
# Run the command.
IJ.run(input_image_plus_copy, "threshold", "")
# Save the ImagePlus object as a new image.
IJ.saveAs(input_image_plus_copy, output_datatype, output_filename)
