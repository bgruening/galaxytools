import sys

from ij import IJ
from ij.plugin import Duplicator

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
input_file = sys.argv[-13]
xleft = int(sys.argv[-12])
width = int(sys.argv[-11])
ytop = int(sys.argv[-10])
height = int(sys.argv[-9])
first_channel = int(sys.argv[-8])
last_channel = int(sys.argv[-7])
first_slice = int(sys.argv[-6])
last_slice = int(sys.argv[-5])
first_frame = int(sys.argv[-4])
last_frame = int(sys.argv[-3])
output_filename = sys.argv[-2]
output_datatype = sys.argv[-1]

# Open the input image file.
input_image_plus = IJ.openImage(input_file)

# Get image dimensions (width, height, nChannels, nSlices, nFrames)
image_dims = input_image_plus.getDimensions()

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()

# Determine if crop in XY is needed:
if xleft != 0 or width != 0 or ytop != 0 or height != 0:
    # Need to define a ROI
    if width == 0:
        width = image_dims[0] - xleft
    if height == 0:
        height = image_dims[1] - ytop
    input_image_plus_copy.setRoi(xleft, ytop, width, height)

# Replace 0's with default:
if last_channel == 0:
    last_channel = image_dims[2]
if last_slice == 0:
    last_slice = image_dims[3]
if last_frame == 0:
    last_frame = image_dims[4]
print("Original dimensions:")
print(image_dims)
# This will also crop in XY is a ROI has been set
input_image_plus_copy = Duplicator().run(input_image_plus_copy, first_channel, last_channel, first_slice, last_slice, first_frame, last_frame)
print("Final dimensions:")
print(input_image_plus_copy.getDimensions())

# Save the ImagePlus object as a new image.
IJ.saveAs(input_image_plus_copy, output_datatype, output_filename)
