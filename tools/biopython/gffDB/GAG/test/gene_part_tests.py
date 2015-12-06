#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.gene_part import GenePart
from src.xrna import XRNA
from src.gene import Gene

class TestGenePart(unittest.TestCase):

    def setUp(self):
        self.gp1 = GenePart()
        self.gp2 = GenePart(feature_type='CDS', indices=[1, 44])
        self.gp2.add_indices([65, 103])
        self.gp2.identifier = ['foo1', 'foo2']
        self.gp2.parent_id = 'mama'
        self.gp3 = GenePart(feature_type='exon')

    def test_constructor(self):
        self.gp1 = GenePart()
        self.assertEquals(0, len(self.gp1.indices))
        self.assertFalse(self.gp1.score)
        self.assertFalse(self.gp1.parent_id)
        self.assertEquals('CDS', self.gp2.feature_type)
        self.assertEquals('+', self.gp1.strand)

    def test_add_indices(self):
        self.assertEquals(2, len(self.gp2.indices))
        self.gp2.add_indices([150, 197])
        self.assertEquals(3, len(self.gp2.indices))

        self.assertEquals(0, len(self.gp3.indices))
        self.gp3.add_indices([77, 144])
        self.assertEquals(1, len(self.gp3.indices))
        # error check
        self.assertRaises(ValueError, self.gp3.add_indices, 7)
        self.assertRaises(ValueError, self.gp3.add_indices, 'foo')

    def test_add_identifier(self):
        self.assertEquals(0, len(self.gp3.identifier))
        self.gp3.add_identifier('7')
        self.assertEquals(1, len(self.gp3.identifier))
    
    def test_add_annotation(self):
        gp = GenePart()
        self.assertFalse(gp.annotations)
        gp.add_annotation("gag_flag", "this gene part rulz")
        self.assertTrue(gp.annotations)

    def test_sort_attributes(self):
        gp = GenePart()
        gp.indices = [[25, 30], [5, 10]] # out of order!
        gp.identifier = ["gp2", "gp1"]
        gp.score = [10, 8]
        self.assertEquals("gp1", gp.identifier[1])
        self.assertEquals([25, 30], gp.indices[0])
        gp.sort_attributes()
        self.assertEquals("gp1", gp.identifier[0])
        self.assertEquals([5, 10], gp.indices[0])

    def test_gagflagged(self):
        gp = GenePart()
        self.assertFalse(gp.gagflagged())
        gp.add_annotation("gag_flag", "awesome flag")
        self.assertTrue(gp.gagflagged())

    def test_length(self):
        self.assertEquals(83, self.gp2.length())
        # what if no indices at all?
        self.assertFalse(self.gp1.length())

    def test_adjust_indices(self):
        self.gp2.adjust_indices(10)
        self.assertEquals(54, self.gp2.indices[0][1])
        # now put it back
        self.gp2.adjust_indices(-10)
        self.assertEquals(65, self.gp2.indices[1][0])

    def test_adjust_indices_after_start_index(self):
        self.gp2.adjust_indices(10, 50)
        self.assertEqual(44, self.gp2.indices[0][1])
        self.assertEqual(75, self.gp2.indices[1][0])

    def test_generate_attribute_entry(self):
        # test .generate_attribute_entry
        expected = "ID=foo1;Parent=mama\n"
        self.assertEquals(expected, self.gp2.generate_attribute_entry(0))
        expected = "ID=foo2;Parent=mama\n"
        self.assertEquals(expected, self.gp2.generate_attribute_entry(1))
        # what if index out of range?
        self.assertFalse(self.gp2.generate_attribute_entry(2))
        # what if no identifier or parent_id?
        self.assertFalse(self.gp1.generate_attribute_entry(0))
        self.gp1.parent_id = 'dad'
        self.assertFalse(self.gp1.generate_attribute_entry(0))

    def test_str(self):
        expected = "CDS (first ID=foo1)"
        self.assertEquals(expected, str(self.gp2))

    def test_to_gff(self):
        # test .to_gff
        seq_name = "sctg_0001_0001"
        source = 'maker'
        strand = '+'
        expected = "sctg_0001_0001\tmaker\tCDS\t1\t44\t.\t+\t.\t"
        expected += "ID=foo1;Parent=mama\n"
        expected += "sctg_0001_0001\tmaker\tCDS\t65\t103\t.\t+\t.\t"
        expected += "ID=foo2;Parent=mama\n"
        actual = self.gp2.to_gff(seq_name=seq_name, source=source)
        self.assertEquals(expected, actual)
        # what if no indices, etc.?
        self.assertFalse(self.gp1.to_gff(seq_name="foo", source="bar"))
        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGenePart))
    return suite

if __name__ == '__main__':
    unittest.main()
