#!/usr/bin/env python3

import argparse
import xml.etree.ElementTree as ET


def local_name(tag):
    return tag.rsplit("}", 1)[-1]


parser = argparse.ArgumentParser()
parser.add_argument("--format", choices=("alto", "page"), required=True)
parser.add_argument("--input", required=True)
parser.add_argument("--output", required=True)
parser.add_argument("--filename", required=True)
args = parser.parse_args()

# Preserve namespace prefixes where possible.
for _, (prefix, uri) in ET.iterparse(args.input, events=("start-ns",)):
    try:
        ET.register_namespace(prefix, uri)
    except ValueError:
        pass

xml_parser = ET.XMLParser(
    target=ET.TreeBuilder(insert_comments=True, insert_pis=True)
)
tree = ET.parse(args.input, parser=xml_parser)
root = tree.getroot()

if args.format == "alto":
    source_information = next(
        element
        for element in root.iter()
        if local_name(element.tag) == "sourceImageInformation"
    )
    filename = next(
        element
        for element in source_information
        if local_name(element.tag) == "fileName"
    )
    filename.text = args.filename

else:
    page = next(
        element
        for element in root.iter()
        if local_name(element.tag) == "Page"
    )
    page.set("imageFilename", args.filename)

tree.write(
    args.output,
    encoding="utf-8",
    xml_declaration=True,
)