import jython_utils
import sys
from ij import IJ
from ij.plugin.filter import Analyzer

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -13 ]
input = sys.argv[ -12 ]
size = sys.argv[ -11 ]
circularity_min = float( sys.argv[ -10 ] )
circularity_max = float( sys.argv[ -9 ] )
show = sys.argv[ -8 ]
display_results = sys.argv[ -7 ]
all_results = sys.argv[ -6 ]
exclude_edges = sys.argv[ -5 ]
include_holes = sys.argv[ -4 ]
tmp_output_path = sys.argv[ -3 ]
output_datatype = sys.argv[ -2 ]
results_path = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )

# Create a copy of the image.
input_image_plus_copy = input_image_plus.duplicate()
image_processor_copy = input_image_plus_copy.getProcessor()
analyzer = Analyzer( input_image_plus_copy )

try:
    # Convert image to binary if necessary.
    if not image_processor_copy.isBinary():
        # Convert the image to binary grayscale.
        IJ.run( input_image_plus_copy, "Make Binary","iterations=1 count=1 edm=Overwrite do=Nothing" )

    # Set the options.
    options = [ 'size=%s' % size ]
    circularity_str = '%.3f-%.3f' % ( circularity_min, circularity_max )
    options.append( 'circularity=%s' % circularity_str )
    if show.find( '_' ) >= 0:
        show_str = '[%s]' % show.replace( '_', ' ' )
    else:
        show_str = show
    options.append( 'show=%s' % show_str ) 
    if jython_utils.asbool( display_results ):
        options.append( 'display' )
        if not jython_utils.asbool( all_results ):
            options.append( 'summarize' )
    if jython_utils.asbool( exclude_edges ):
        options.append( 'exclude' )
    if jython_utils.asbool( include_holes ):
        options.append( 'include' )
    # Always run "in_situ".
    options.append( 'in_situ' )

    # Run the command.
    IJ.run( input_image_plus_copy, "Analyze Particles...", " ".join( options ) )
    
    # Save outputs.
    if tmp_output_path not in [ None, 'None' ]:
        # Save the ImagePlus object as a new image.
        IJ.saveAs( input_image_plus_copy, output_datatype, tmp_output_path )
    if jython_utils.asbool( display_results ) and results_path not in [ None, 'None' ]:
        results_table = analyzer.getResultsTable()
        results_table.saveAs( results_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
