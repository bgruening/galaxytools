import sys

from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
input = sys.argv[-7]
iterations = int(sys.argv[-6])
count = int(sys.argv[-5])
black_background = sys.argv[-4] == "yes"
pad_edges_when_eroding = sys.argv[-3] == "yes"
tmp_output_path = sys.argv[-2]
output_datatype = sys.argv[-1]

# Open the input image file.
input_image_plus = IJ.openImage(input)

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()

# Set binary options.
options = ["edm=Overwrite", "iterations=%d" % iterations, "count=%d" % count]
if pad_edges_when_eroding:
    options.append("pad")
if black_background:
    options.append("black")
options = " ".join(options)
IJ.run(input_image_plus_copy, "Options...", options)

# Run the command.
IJ.run(input_image_plus_copy, "Make Binary", "")

# Save the ImagePlus object as a new image.
IJ.saveAs(input_image_plus_copy, output_datatype, tmp_output_path)
