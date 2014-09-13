#!/usr/bin/env python

"""
    Create automatic tool_conf.xml file.
    Add this file to tool_config_file in your universe_wsgi.conf.
    For example:
    tool_config_file = tool_conf.xml,shed_tool_conf.xml,chemicaltoolbox.xml
"""

import os
import sys
import argparse


def main(options):
    ohandle = options.outfile
    ohandle.write("""<?xml version="1.0"?>\n<toolbox tool_path="%s">\n""" % options.rel_tool_path)
    for root, dirs, files in os.walk( options.input_path, topdown=True):
        for dirname in dirs:
            if dirname.startswith('.'):
                dirs.remove( dirname )

        tools = list()
        for filename in files:
            if filename in ['tool_dependencies.xml', 'tool_conf.xml', 'datatypes_conf.xml']:
                continue
            path = os.path.join( root, filename )
            if os.path.splitext( path )[-1] == '.xml':
                tools.append('\t\t<tool file="%s"/>\n' % os.path.abspath( os.path.join(root,filename) ))
        if tools:
            ohandle.write('\t<section id="%s" name="%s" version="">\n' % ('bgruening_%s' % os.path.basename(root), os.path.basename(root)))
            ohandle.write(''.join(tools))
            ohandle.write('\t</section>\n')
    ohandle.write('</toolbox>\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create tool_conf.xml file.')

    parser.add_argument("-i", "--input", dest="input_path",
                    required=True,
                    help="input directory with to parse recursivly")
    parser.add_argument("-t", dest="rel_tool_path",
                    help="relative tool path", required=True)
    parser.add_argument('-o', '--outfile', type=argparse.FileType('w'),
                     default=sys.stdout)

    options = parser.parse_args()

    main(options)
