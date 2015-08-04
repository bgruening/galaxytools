import jython_utils
import math
import sys
from ij import IJ
from ij import ImagePlus
from skeleton_analysis import AnalyzeSkeleton_

BASIC_NAMES = [ 'Branches', 'Junctions', 'End-point Voxels',
                'Junction Voxels', 'Slab Voxels', 'Average branch length',
                'Triple Points', 'Quadruple Points', 'Maximum Branch Length' ]
DETAIL_NAMES = [ 'Skeleton ID', 'Branch length', 'V1 x', 'V1 y', 'V1 z', 'V2 x',
                 'V2 y', 'V2 z', 'Euclidean distance' ]
VALID_IMAGE_TYPES = [ ImagePlus.GRAY8 ]

def asbool( val ):
    return val == 'yes'

def get_euclidean_distance( vertex1, vertex2 ):
    x1, y1, z1 = get_points( vertex1 )
    x2, y2, z2 = get_points( vertex2 )
    return math.sqrt( math.pow( ( x2 - x1 ), 2 ) +
                      math.pow( ( y2 - y1 ), 2 ) +
                      math.pow( ( z2 - z1 ), 2 ) )

def get_graph_length( graph ):
    length = 0
    for edge in graph.getEdges():
        length = length + edge.getLength()
    return length

def get_points( vertex ):
    # An array of Point, which has attributes x,y,z.
    point = vertex.getPoints()[ 0 ]
    return point.x, point.y, point.z
    
def get_sorted_edge_lengths( graph ):
    # Return graph edges sorted from longest to shortest.
    edges = graph.getEdges()
    edges = sorted( edges, key=lambda edge: edge.getLength(), reverse=True )
    return edges

def get_sorted_graph_lengths( result ):
    # Get the separate graphs (skeletons).
    graphs = result.getGraph()
    # Sort graphs from longest to shortest.
    graphs = sorted( graphs, key=lambda g: get_graph_length( g ), reverse=True )
    return graphs

def save( result, output, show_detailed_info, calculate_largest_shortest_path, sep='\t' ):
    num_trees = int( result.getNumOfTrees() )
    outf = open( output, 'wb' )
    outf.write( '# %s\n' % sep.join( BASIC_NAMES ) )
    for index in range( num_trees ):
        outf.write( '%d%s' % ( result.getBranches()[ index ], sep ) )
        outf.write( '%d%s' % ( result.getJunctions()[ index ], sep ) )
        outf.write( '%d%s' % ( result.getEndPoints()[ index ], sep ) )
        outf.write( '%d%s' % ( result.getJunctionVoxels()[ index ], sep ) )
        outf.write( '%d%s' % ( len( result.getListOfSlabVoxels() ), sep ) )
        outf.write( '%.3f%s' % ( result.getAverageBranchLength()[ index ], sep ) )
        outf.write( '%d%s' % ( result.getTriples()[ index ], sep ) )
        outf.write( '%d%s' % ( result.getQuadruples()[ index ], sep ) )
        outf.write( '%.3f' % result.getMaximumBranchLength()[ index ] )
        if calculate_largest_shortest_path:
            outf.write( '%s%.3f%s' % ( sep, result.shortestPathList.get( index ), sep ) )
            outf.write( '%d%s' % ( result.spStartPosition[ index ][ 0 ], sep ) )
            outf.write( '%d%s' % ( result.spStartPosition[ index ][ 1 ], sep ) )
            outf.write( '%d\n' % result.spStartPosition[ index ][ 2 ] )
        else:
            outf.write( '\n' )
    if show_detailed_info:
        outf.write( '# %s\n' % sep.join( DETAIL_NAMES ) )
        # The following index is a placeholder for the skeleton ID.
        # The terms "graph" and "skeleton" refer to the same thing.
        # Also, the SkeletonResult.java code states that the
        # private Graph[] graph object is an array of graphs (one
        # per tree).
        for index, graph in enumerate( get_sorted_graph_lengths( result ) ):
            for edge in get_sorted_edge_lengths( graph ):
                vertex1 = edge.getV1()
                x1, y1, z1 = get_points( vertex1 )
                vertex2 = edge.getV2()
                x2, y2, z2 = get_points( vertex2 )
                outf.write( '%d%s' % ( index+1, sep ) )
                outf.write( '%.3f%s' % ( edge.getLength(), sep ) )
                outf.write( '%d%s' % ( x1, sep ) )
                outf.write( '%d%s' % ( y1, sep ) )
                outf.write( '%d%s' % ( z1, sep ) )
                outf.write( '%d%s' % ( x2, sep ) )
                outf.write( '%d%s' % ( y2, sep ) )
                outf.write( '%d%s' % ( z2, sep ) )
                outf.write( '%.3f' % get_euclidean_distance( vertex1, vertex2 ) )
                if calculate_largest_shortest_path:
                    # Keep number of separated items the same for each line.
                    outf.write( '%s %s' % ( sep, sep ) )
                    outf.write( ' %s' % sep )
                    outf.write( ' %s' % sep )
                    outf.write( ' \n' )
                else:
                    outf.write( '\n' )
    outf.close()

# Fiji Jython interpreter implements Python 2.5 which does not
# provide support for argparse.
error_log = sys.argv[ -7 ]
input = sys.argv[ -6 ]
prune_cycle_method = sys.argv[ -5 ]
prune_ends = asbool( sys.argv[ -4 ] )
calculate_largest_shortest_path = asbool( sys.argv[ -3 ] )
if calculate_largest_shortest_path:
    BASIC_NAMES.extend( [ 'Longest Shortest Path', 'spx', 'spy', 'spz' ] )
    DETAIL_NAMES.extend( [ ' ', ' ', ' ', ' ' ] )
show_detailed_info = asbool( sys.argv[ -2 ] )
output = sys.argv[ -1 ]

# Open the input image file.
input_image_plus = IJ.openImage( input )
bit_depth = input_image_plus.getBitDepth()
image_type = input_image_plus.getType()

# Create a copy of the image.
input_image_plus_copy = input_image_plus.createImagePlus()
image_processor_copy = input_image_plus.getProcessor().duplicate()
input_image_plus_copy.setProcessor( "iCopy", image_processor_copy )

try:
    if image_type not in VALID_IMAGE_TYPES:
        # Convert the image to binary grayscale.
        IJ.run( input_image_plus_copy, "Make Binary","iterations=1 count=1 edm=Overwrite do=Nothing" )
    # Run AnalyzeSkeleton
    analyze_skeleton = AnalyzeSkeleton_()
    analyze_skeleton.setup( "", input_image_plus_copy )
    if prune_cycle_method == 'none':
        prune_index  = analyze_skeleton.NONE
    elif prune_cycle_method == 'shortest_branch':
        prune_index  = analyze_skeleton.SHORTEST_BRANCH
    elif prune_cycle_method == 'lowest_intensity_voxel':
        prune_index  = analyze_skeleton.LOWEST_INTENSITY_VOXEL
    elif prune_cycle_method == 'lowest_intensity_branch':
        prune_index  = analyze_skeleton.LOWEST_INTENSITY_BRANCH
    result = analyze_skeleton.run( prune_index,
                                   prune_ends,
                                   calculate_largest_shortest_path,
                                   input_image_plus_copy,
                                   True,
                                   True )
    # Save the results.
    save( result, output, show_detailed_info, calculate_largest_shortest_path )
except Exception, e:
    jython_utils.handle_error( error_log, str( e ) )
