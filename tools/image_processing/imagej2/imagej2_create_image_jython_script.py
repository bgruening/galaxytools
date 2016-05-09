import sys
from ij import IJ

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
title = sys.argv[ -6 ]
width = int( sys.argv[ -5 ] )
height = int( sys.argv[ -4 ] )
depth = int( sys.argv[ -3 ] )
type = sys.argv[ -2 ].replace( '_', ' ' )
tmp_image_path = sys.argv[ -1 ]

imp = IJ.newImage( title, type, width, height, depth )
IJ.save( imp, "%s" % tmp_image_path )
