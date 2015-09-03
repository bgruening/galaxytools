Galaxy wrappers for ImageJ2 tools
==================================

ImageJ2 is a new version of ImageJ for the next generation of multidimensional image data, with a focus on scientific imaging. Its central goal is to broaden the paradigm of ImageJ beyond the limitations of ImageJ 1.x, to support the next generation of multidimensional scientific imaging.

Fiji is an image processing package. It can be described as a "batteries-included" distribution of ImageJ (and ImageJ2), bundling Java, Java3D and a lot of plugins organized into a coherent menu structure. Fiji compares to ImageJ as Ubuntu compares to Linux.

More informations is available at:

* [http://fiji.sc/ImageJ2](http://fiji.sc/ImageJ2)
* [http://fiji.sc/Fiji](http://fiji.sc/Fiji)


Installation
============

Galaxy tool wrappers use specified Fiji Lifeline versions available from [http://fiji.sc/Downloads](http://fiji.sc/Downloads).  Galaxy should be able to automatically install this package.

The wrappers are available at [https://github.com/bgruening/galaxytools/tree/master/tools/image_processing/imagej2](https://github.com/bgruening/galaxytools/tree/master/tools/image_processing/imagej2).


Use Docker
==========

A docker image that installs Galaxy with these imaging tools is available at [https://github.com/bgruening/galaxy-imaging](https://github.com/bgruening/galaxy-imaging).


Using Fiji with Galaxy tools
============================

Galaxy ImageJ2 tool wrappers generate a command line that calls a Python script, passing it a series of arguments including a Jython script named jython_script.py that resides in the same directory as the tool wrapper.  During tool execution, the Python script will call ImageJ2 with the --headless argument to run without the ImageJ2 GUI.  The Jython script is also passed to ImageJ2 along with all command line arguments that it expects.  ImageJ2 will execute the Jython script, passing the expected arguments.  The command line to run ImageJ2 from a Galaxy tool wrapper looks something like this:

`ImageJ2 --ij2 --headless --jython ~jython_script.py arg1, arg2, ...`

Each tool execution starts the ImageJ2 application within a Java virtual machine (JVM).  When ImageJ2 is finished processing the Jython script, the results are either written to a file or returned to the calling Galaxy process.  The JVM is shut down, and the Galaxy job terminates.  This approach provides the ability to run ImageJ2 tools from Galaxy on any supported HPC environment.

Of course, eliminating the ImageJ2 GUI restricts us to wrapping only those ImageJ2 plugins that do not require any GUI components (i.e., the ImageJ2 window manager).  Plugins are written by an open community, so not all of them are written in such a way that they can be executed from the command line and produce useful results.  For example, some plugins create one or more images that can only be accessed via calls to the ImageJ2 window manager, and running in headless mode eliminates the window manager as well as other GUI components.

Those familiar with ImageJ2 will find differences with this general pattern for executing ImageJ2 tools within Galaxy.  ImageJ2 accounts for user defined global preferences which are available to tools throughout the session, and an image can be uploaded and run through any number of available tools, saving only the final image.  While Galaxy currently does not account for user preferences defined in ImageJ2, enhancements to the Galaxy framework are planned that will accomodate these kinds of settings (e.g., binary image options).  Also, since Galaxy initiates a new ImageJ2 session with each tool execution, initial images are uploaded to ImageJ2 and resulting images are saved for each tool execution.

The Galaxy ImageJ2 tools currently fall into the following categories.  Additional tools will be added at a steady pace.

Working with Pixels
===================
These Galaxy tools wrap the ImageJ2 plugins that are available in the ImageJ2 Process → Math menu.

* **Operate on pixels** - Applies a mathematical expression (add, subtract, multiply, etc.) to each pixel in the image.  When the resulting pixel value overflows/underflows the legal range of the image's data type, the value is reset to the maximum/minimum value.

Filters and Effects
===================
These Galaxy tools wrap the ImageJ2 plugins that are available in the ImageJ2 Process menu.

* **Smooth** - Blurs the image by replacing each pixel with the average of its 3x3 neighborhood.
* **Sharpen** - Increases contrast and accentuates detail in the image, but may also accentuate noise.
* **Find Edges** - Uses a Sobel edge detector to highlight sharp changes in intensity in the active image.
* **Add shadow effect** - Produces a shadow effect, with light appearing to come from the selected direction (East, North, Northeast, Northwest, South, Southeast, Southwest and West).
* **Find Maxima** - Determines the local maxima in an image and creates a binary (mask-like) image of the same size with the maxima (or one segmented particle per maximum) marked.
* **Enhance contrast** - Enhances image contrast by using either normalization (contrast stretching) or histogram equalization.
* **Add or remove noise** - Adds specified noise to or removes noise from images.

Binary Image Tools
==================
These Galaxy tools wrap the ImageJ2 plugins that are available in the ImageJ2 Process → Binary menu.

* **Convert to binary** - Converts an image into a binary (black and white) image.
* **Adjust threshold** - Sets lower and upper threshold values, segmenting grayscale images into
features of interest and background.
* **Watershed segmentation** - Automatically separates or cuts apart particles that touch.
* **Analyze particles** - Analyzes the particles in a segmented binary image, providing information about
each particle in the image.
* **Skeletonize images** - Uses the Skeletonize3D plugin to find the centerlines (”skeleton”) of objects in the input image.  Skeletonize3d is a plugin written by Ignacio Arganda-Carreras that offers several advantages over the legacy skeletonization algorithm of ImageJ available in the Process → Binary → Skeletonize menu item.  Skeletonize works only with binary 2D images.  Skeletonize3D works with 8-bit 2D images and stacks, expecting the image to be binary. If not, Skeletonize3D considers all pixel values above 0 to be white (255).  While Skeletonize↑ relies on Black background value, the output of Skeletonize3D always has a value of 255 at the skeleton and 0 at background pixels, independently of the Black background option.
* **Analyze skeleton** - Tags all pixel/voxels in a skeleton image and then counts all its junctions,
triple and quadruple points and branches, and measures their average and maximum length.
* **Convert binary image to EDM** - Converts a binary image into a 8-bit grayscale Euclidean Distance Map (EDM). Each foreground (nonzero) pixel in the binary image is assigned a value equal to its distance from the nearest background (zero) pixel.

**Interpreting binary Images in ImageJ2**

Binary images are thresholded to only two values, typically 0 and 1, but often — as with ImageJ — 0 and 255, that represent black and white on an 8-bit scale.

The interpretation of binary images is not universal. While some software packages will always perform binary operations on 255 values (or 1, or any non-zero value), ImageJ takes into account the foreground and background colors of the binary image.

In ImageJ, the **Black background** global preference setting defines not only how new binary images will be created, but also how previously created images are interpreted. This means objects will be inferred on a image-per-image basis.  As such, inverting the LUT (i.e., pixels with a value of zero are white and pixels with a value 255 are black) of a binary image without updating the black background option may lead to unexpected results.  This issue can currently be avoided by properly selecting the **Black background** option available on all Galaxy binary image tools.

BunwarpJ Plugin Tools
=====================
These Galaxy tools wrap the bUnwarpJ plugin [http://fiji.sc/BUnwarpJ](http://fiji.sc/BUnwarpJ).

* **Adapt an elastic transformation** - Adapts an elastic transformation to a new image size by transforming the
coefficients of a specified elastic transformation according to a real image factor.
* **Align two images** - Performs a simultaneous registration of two images, A and B. Image A is elastically deformed
in order to look as similar as possible to image B, and, at the same time, the "inverse"
transformation (from B to A) is also calculated so a pseudo-invertibility of the final deformation
could be guaranteed.  Two images are produced: the deformed versions of A and B images.
* **Compare opposite elastic deformations** - Calculates the warping index of two opposite elast transformations, i.e. the average of the geometrical distance between every pixel and its version after applying both transformations (direct and inverse).
* **Compare elastic and raw deformation** - Calculates the warping index of an elastic transformation and a raw transformation.
* **Compare two raw deformations** - Calculates the warping index of two raw transformations (same direction).
* **Compose two elastic transformations** - Composes two elastic transformations into a raw transformation.
* **Compose two raw transformations** - Composes two raw transformations into another raw transformation.
* **Compose a raw and an elastic transformation** - Composes a raw transformation and an elastic transformation
into a raw transformation.
* **Convert elastic transformation to raw** - Converts an elastic (i.e., B-spline ) transformation file into a raw transformation file.
* **Apply elastic transformation** - Applies an elastic transformation to an image, producing another image which is elastically
deformed according to the transformation.
* **Apply raw transformation** - Applies a raw transformation to an image, producing another image which is deformed according
to the transformation.

Other Tools
===========
* **Create new image** - Creates a new image of a selected type, size, depth and format.
* **Convert image format** - Converts the format of an input image file, producing an output image.

Licence
=======

Fiji is released as open source under the GNU General Public License: [http://www.gnu.org/licenses/gpl.html](http://www.gnu.org/licenses/gpl.html)

Fiji builds on top of the ImageJ2 core, which is licensed under the permissive BSD 2-Clause license: [http://opensource.org/licenses/BSD-2-Clause](http://opensource.org/licenses/BSD-2-Clause)

Plugins and other components have their own licenses.

