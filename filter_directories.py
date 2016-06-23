#!/usr/bin/env python

import os
import sys
import lxml.etree as etree

def is_tool(path):
    xml = etree.parse(path)
    return xml.getroot().tag == 'tool'


with open('.tt_blacklist') as handle:
    bl = [tool.strip() for tool in handle]

for directory in sys.stdin:
    directory = directory.strip()
    if directory.startswith('packages') or directory.startswith('data_managers'):
        continue
    while directory:
        if os.path.exists( os.path.join(directory, '.shed.yml') ) or os.path.exists( os.path.join(directory, '.shed.yaml') ):
            if directory not in bl:
                for filename in os.listdir(directory):
                    xml_path = os.path.join(directory, filename)
                    if xml_path.endswith('.xml') and is_tool(xml_path):
                        print(xml_path)
            break
        else:
            directory = os.path.dirname( directory )
