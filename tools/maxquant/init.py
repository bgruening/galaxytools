#!/usr/bin/env python3
"""Initialize MaxQuant tool for use with a new version of
modifications/enzymes.xml.

TODO: Append function: only add modifications that are not
already present, add modification entries to conda maxquant

Authors: Damian Glaetzer <d.glaetzer@mailbox.org>
"""

import xml.etree.ElementTree as ET
import os
import sys
import re
from xml.dom import minidom

usage = '\n'.join(("Usage: {} MODS_FILE ENZYMES_FILE".format(sys.argv[0]),
                   "FILES are the modifications/enzymes.xml of MaxQuant.",
                   "Updates modification parameters in macros.xml."))


if len(sys.argv) != 3 or not (os.path.isfile(sys.argv[1])
                              and os.path.isfile(sys.argv[2])):
    print(usage)

mods_root = ET.parse(sys.argv[1]).getroot()

mods = mods_root.findall('modification')
standard_mods = []
label_mods = []
for m in mods:
    if m.findtext('type') == 'Standard':
        standard_mods.append(m.get('title'))
    elif m.findtext('type') == 'Label':
        label_mods.append(m.get('title'))


enzymes_root = ET.parse(sys.argv[2]).getroot()

enzymes = enzymes_root.findall('enzyme')
enzymes_list = [e.get('title') for e in enzymes]

macros_root = ET.parse('./macros.xml').getroot()
for child in macros_root:
    if child.get('name') == 'modification':
        child.clear()
        child.tag = 'xml'
        child.set('name', 'modification')
        for m in standard_mods:
            ET.SubElement(child, 'expand', attrib={'macro': 'mod_option',
                                                   'value': m})
    elif child.get('name') == 'label':
        child.clear()
        child.tag = 'xml'
        child.set('name', 'label')
        for m in label_mods:
            ET.SubElement(child, 'expand', attrib={'macro': 'mod_option',
                                                   'value': m})

    elif child.get('name') == 'proteases':
        child.clear()
        child.tag = 'xml'
        child.set('name', 'proteases')
        for e in enzymes_list:
            ET.SubElement(child, 'expand', attrib={'macro': 'mod_option',
                                                   'value': e})

rough_string = ET.tostring(macros_root, 'utf-8')
reparsed = minidom.parseString(rough_string)
pretty = reparsed.toprettyxml(indent="\t")
even_prettier = re.sub(r"\n\s+\n", r"\n", pretty)
with open('./macros.xml', 'w') as f:
    print(even_prettier, file=f)
