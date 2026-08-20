#!/usr/bin/env python

import unittest
from mock import Mock
from src.seq_helper import SeqHelper

class TestSeqHelper(unittest.TestCase):

    def setUp(self):
        self.helper = SeqHelper("nnnGATTACAnnnnnNNNNNgattacaNNN")

    def test_initialize(self):
        self.assertEquals("nnnGATTACAnnnnnNNNNNgattacaNNN", self.helper.full_sequence)

    def test_mrna_to_fasta(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        mrna.exon = Mock()
        mrna.exon.indices = [[4, 10], [21, 27]]
        mrna.strand = '+'
        expected = ">foo_mrna\nGATTACAgattaca\n"
        self.assertEquals(expected, self.helper.mrna_to_fasta(mrna))

    def test_mrna_to_fasta_reverse_strand(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        mrna.exon = Mock()
        mrna.exon.indices = [[4, 10], [21, 27]]
        mrna.strand = '-'
        expected = ">foo_mrna\ntgtaatcTGTAATC\n"
        self.assertEquals(expected, self.helper.mrna_to_fasta(mrna))

    def test_mrna_to_cds_fasta(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        mrna.cds = Mock()
        mrna.cds.indices = [[4, 10], [21, 27]]
        mrna.strand = '-'
        expected = ">foo_mrna CDS\ntgtaatcTGTAATC\n"
        self.assertEquals(expected, self.helper.mrna_to_cds_fasta(mrna))

    def test_mrna_to_protein_fasta(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        mrna.cds = Mock()
        mrna.cds.indices = [[4, 10], [21, 27]]
        mrna.strand = '+'
        expected = ">foo_mrna protein\nDYRL\n"
        self.assertEquals(expected, self.helper.mrna_to_protein_fasta(mrna))

    def test_mrna_to_protein_fasta_reverse(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        mrna.cds = Mock()
        mrna.cds.indices = [[21, 27], [4, 10]]
        mrna.strand = '-'
        expected = ">foo_mrna protein\nCNL*\n"
        self.assertEquals(expected, self.helper.mrna_to_protein_fasta(mrna))

    def test_mrna_contains_internal_stop(self):
        helper = SeqHelper("gattacaTAGgattaca") # TAG = stop codon
        mrna = Mock()
        mrna.strand = '+'
        mrna.cds = Mock()
        mrna.cds.indices = [[2, 4], [8, 14]]
        # verify that the internal stop is in the sequence ...
        self.assertEquals("attTAGgatt", helper.get_sequence_from_indices('+', mrna.cds.indices))
        # make sure the method catches it
        self.assertTrue(helper.mrna_contains_internal_stop(mrna))

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSeqHelper))
    return suite

if __name__ == '__main__':
    unittest.main()
