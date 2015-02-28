#!/usr/bin/env python

import sys
import re
import tempfile
#from rpy import *
import rpy2.robjects as robjects
r = robjects.r
from rpy2.robjects.vectors import DataFrame

assert sys.version_info[:2] >= ( 2, 4 )

def stop_err( msg ):
    sys.stderr.write( msg )
    sys.exit()

def S3_METHODS( all="key" ):
    Group_Math =  [ "abs", "sign", "sqrt", "floor", "ceiling", "trunc", "round", "signif",
        "exp", "log", "cos", "sin", "tan", "acos", "asin", "atan", "cosh", "sinh", "tanh",
        "acosh", "asinh", "atanh", "lgamma", "gamma", "gammaCody", "digamma", "trigamma",
        "cumsum", "cumprod", "cummax", "cummin", "c" ]
    Group_Ops = [ "+", "-", "*", "/", "^", "%%", "%/%", "&", "|", "!", "==", "!=", "<", "<=", ">=", ">", "(", ")", "~", "," ]
    if all is "key":
        return { 'Math' : Group_Math, 'Ops' : Group_Ops }

def main():
    try:
        datafile = sys.argv[1]
        outfile_name = sys.argv[2]
        expression = sys.argv[3]
    except: 
        stop_err( 'Usage: python gsummary.py input_file ouput_file expression' )

    math_allowed = S3_METHODS()[ 'Math' ]
    ops_allowed = S3_METHODS()[ 'Ops' ]

    # Check for invalid expressions
    for word in re.compile( '[a-zA-Z]+' ).findall( expression ):
        if word and not word in math_allowed: 
            stop_err( "Invalid expression '%s': term '%s' is not recognized or allowed" %( expression, word ) )
    symbols = set()
    for symbol in re.compile( '[^a-z0-9\s]+' ).findall( expression ):
        if symbol and not symbol in ops_allowed:
            stop_err( "Invalid expression '%s': operator '%s' is not recognized or allowed" % ( expression, symbol ) )
        else:
            symbols.add( symbol )
    if len( symbols ) == 1 and ',' in symbols:
        # User may have entered a comma-separated list r_data_frame columns
        stop_err( "Invalid columns '%s': this tool requires a single column or expression" % expression )

    # Find all column references in the expression
    cols = []
    for col in re.compile( 'c[0-9]+' ).findall( expression ):
        try:
            cols.append( int( col[1:] ) - 1 )
        except:
            pass
 
    tmp_file = tempfile.NamedTemporaryFile( 'w+b' )
    # Write the R header row to the temporary file
    hdr_str = "\t".join( "c%s" % str( col+1 ) for col in cols )
    tmp_file.write( "%s\n" % hdr_str )
    skipped_lines = 0
    first_invalid_line = 0
    i = 0
    for i, line in enumerate( file( datafile ) ):
        line = line.rstrip( '\r\n' )
        if line and not line.startswith( '#' ):
            valid = True
            fields = line.split( '\t' )
            # Write the R data row to the temporary file
            for col in cols:
                try:
                    float( fields[ col ] )
                except:
                    skipped_lines += 1
                    if not first_invalid_line:
                        first_invalid_line = i + 1
                    valid = False
                    break
            if valid:
                data_str = "\t".join( fields[ col ] for col in cols )
                tmp_file.write( "%s\n" % data_str )
    tmp_file.flush()

    if skipped_lines == i + 1:
        stop_err( "Invalid column or column data values invalid for computation.  See tool tips and syntax for data requirements." )
    else:
        # summary function and return labels
        summary_func = r( "function( x ) { c( sum=sum( as.numeric( x ), na.rm=T ), mean=mean( as.numeric( x ), na.rm=T ), stdev=sd( as.numeric( x ), na.rm=T ), quantile( as.numeric( x ), na.rm=TRUE ) ) }" )
        headings = [ 'sum', 'mean', 'stdev', '0%', '25%', '50%', '75%', '100%' ]
        headings_str = "\t".join( headings )
        
        #r.set_default_mode( NO_CONVERSION )
        #r_data_frame = r.read_table( tmp_file.name, header=True, sep="\t" )
        r_data_frame = DataFrame.from_csvfile( tmp_file.name, header=True, sep="\t" )
        
        outfile = open( outfile_name, 'w' )

        for col in re.compile( 'c[0-9]+' ).findall( expression ):
            r.assign( col, r[ "$" ]( r_data_frame, col ) )
        try:
            summary = summary_func( r( expression ) )
        except RException, s:
            outfile.close()
            stop_err( "Computation resulted in the following error: %s" % str( s ) )
        #summary = summary.as_py( BASIC_CONVERSION )
        outfile.write( "#%s\n" % headings_str )
        print summary
        print summary.r_repr()
        outfile.write( "%s\n" % "\t".join( [ "%g" % ( summary.rx2( k )[0] ) for k in headings ] ) )
        outfile.close()

        if skipped_lines:
            print "Skipped %d invalid lines beginning with line #%d.  See tool tips for data requirements." % ( skipped_lines, first_invalid_line )        

if __name__ == "__main__": main()
