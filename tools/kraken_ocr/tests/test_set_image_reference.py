from pathlib import Path
from subprocess import run
from sys import executable
from tempfile import TemporaryDirectory
from unittest import TestCase

from lxml import etree


class ImageReferenceTest(TestCase):
    def test_preserves_interleaved_reading_order(self):
        with TemporaryDirectory() as directory:
            tmp_path = Path(directory)
            source = tmp_path / "source.xml"
            source.write_text(
                """<alto xmlns="http://www.loc.gov/standards/alto/ns-v4#">
                <Description><sourceImageInformation><fileName>old.png</fileName></sourceImageInformation></Description>
                <ReadingOrder><OrderedGroup ID="outer">
                    <ElementRef REF="A"/>
                    <OrderedGroup ID="inner"><ElementRef REF="B"/></OrderedGroup>
                    <ElementRef REF="C"/>
                </OrderedGroup></ReadingOrder>
            </alto>"""
            )
            destination = tmp_path / "output.xml"
            script = Path(__file__).resolve().parents[1] / "set_image_reference.py"
            run([executable, script, source, "input.jpg", destination], check=True)
            tree = etree.parse(destination)
            self.assertEqual(tree.findtext(".//{*}fileName"), "input.jpg")
            self.assertEqual(tree.xpath("//@REF"), ["A", "B", "C"])
