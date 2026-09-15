"""Point a staged ALTO/PAGE annotation at its paired Galaxy image."""
import sys
import xml.etree.ElementTree as ET


def stage_annotation(source, destination, image_name):
    tree = ET.parse(source)
    for element in tree.iter():
        local_name = element.tag.rsplit('}', 1)[-1]
        if local_name == 'fileName':
            element.text = image_name
        elif local_name == 'Page' and 'imageFilename' in element.attrib:
            element.set('imageFilename', image_name)
    tree.write(destination, encoding='utf-8', xml_declaration=True)


if __name__ == '__main__':
    stage_annotation(*sys.argv[1:])
