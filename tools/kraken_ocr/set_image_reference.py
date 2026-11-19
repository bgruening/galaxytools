import sys

from lxml import etree


def set_image_reference(source, image_filename, destination):
    tree = etree.parse(source, etree.XMLParser(resolve_entities=False))
    root = tree.getroot()
    name = etree.QName(root)
    ns = f"{{{name.namespace}}}" if name.namespace else ""
    if name.localname == "alto":
        reference = root.find(f"{ns}Description/{ns}sourceImageInformation/{ns}fileName")
        if reference is None:
            raise ValueError("Missing image reference in segmentation XML")
        reference.text = image_filename
    elif name.localname == "PcGts":
        page = root.find(f"{ns}Page")
        if page is None or "imageFilename" not in page.attrib:
            raise ValueError("Missing image reference in segmentation XML")
        page.set("imageFilename", image_filename)
    else:
        raise ValueError("Expected an ALTO or PAGE segmentation XML")
    tree.write(destination, encoding="UTF-8", xml_declaration=True)


if __name__ == "__main__":
    set_image_reference(*sys.argv[1:])
