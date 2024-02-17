import sys

from ij import IJ
from ij.plugin.filter import Analyzer

OPTIONS = ["edm=Overwrite", "iterations=1", "count=1"]

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[-14]
input_file = sys.argv[-13]
black_background = sys.argv[-12] == "yes"
size = sys.argv[-11]
circularity_min = float(sys.argv[-10])
circularity_max = float(sys.argv[-9])
show = sys.argv[-8]
display_results = sys.argv[-7] == "yes"
all_results = sys.argv[-6] == "yes"
exclude_edges = sys.argv[-5] == "yes"
include_holes = sys.argv[-4] == "yes"
output_filename = sys.argv[-3]
output_datatype = sys.argv[-2]
results_path = sys.argv[-1]

# Open the input image file.
input_image_plus = IJ.openImage(input_file)

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()
analyzer = Analyzer(input_image_plus_copy)

# Set binary options.
options_list = OPTIONS
if black_background:
    options_list.append("black")
options = " ".join(options_list)
IJ.run(input_image_plus_copy, "Options...", options)

if not image_processor_copy.isBinary():
    # Convert the image to binary grayscale.
    IJ.run(input_image_plus_copy, "Make Binary", "")

# Set the options.
options = ["size=%s" % size]
circularity_str = "%.3f-%.3f" % (circularity_min, circularity_max)
options.append("circularity=%s" % circularity_str)
if show.find("_") >= 0:
    show_str = "[%s]" % show.replace("_", " ")
else:
    show_str = show
options.append("show=%s" % show_str)
if display_results:
    options.append("display")
    if not all_results:
        options.append("summarize")
if exclude_edges:
    options.append("exclude")
if include_holes:
    options.append("include")
# Always run "in_situ".
options.append("in_situ")

# Run the command.
IJ.run(input_image_plus_copy, "Analyze Particles...", " ".join(options))

# Save outputs.
if len(output_filename) > 0:
    # Save the ImagePlus object as a new image.
    IJ.saveAs(input_image_plus_copy, output_datatype, output_filename)
if display_results and len(results_path) > 0:
    results_table = analyzer.getResultsTable()
    results_table.saveAs(results_path)
