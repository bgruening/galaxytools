import sys

from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[-8]
input_file = sys.argv[-7]
iterations = int(sys.argv[-6])
count = int(sys.argv[-5])
black_background = sys.argv[-4] == "yes"
pad_edges_when_eroding = sys.argv[-3] == "yes"
output_filename = sys.argv[-2]
output_datatype = sys.argv[-1]

# Open the input image file.
input_image_plus = IJ.openImage(input_file)

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()

# Set binary options.
options_list = ["edm=Overwrite", "iterations=%d" % iterations, "count=%d" % count]
if black_background:
    options_list.append("black")
if pad_edges_when_eroding:
    options_list.append("pad")
options = " ".join(options_list)
IJ.run(input_image_plus_copy, "Options...", options)

# Convert image to binary if necessary.
if not image_processor_copy.isBinary():
    # Convert the image to binary grayscale.
    IJ.run(input_image_plus_copy, "Make Binary", "")

# Run the command.
IJ.run(input_image_plus_copy, "Distance Map", "")
# Save the ImagePlus object as a new image.
IJ.saveAs(input_image_plus_copy, output_datatype, output_filename)
