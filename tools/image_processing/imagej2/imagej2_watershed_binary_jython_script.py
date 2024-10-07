import sys

from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
input = sys.argv[-4]
black_background = sys.argv[-3] == "yes"
tmp_output_path = sys.argv[-2]
output_datatype = sys.argv[-1]

# Open the input image file.
input_image_plus = IJ.openImage(input)

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()

# Set binary options.
options = ["edm=Overwrite", "iterations=1", "count=1"]
if black_background:
    options.append("black")
options = " ".join(options)
IJ.run(input_image_plus_copy, "Options...", options)

# Convert image to binary if necessary.
if not image_processor_copy.isBinary():
    # Convert the image to binary grayscale.
    IJ.run(input_image_plus_copy, "Make Binary", "")

# Run the command.
IJ.run(input_image_plus_copy, "Watershed", "stack")

# Save the ImagePlus object as a new image.
IJ.saveAs(input_image_plus_copy, output_datatype, tmp_output_path)
