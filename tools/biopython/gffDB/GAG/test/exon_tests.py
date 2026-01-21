#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.gene_part import GenePart, get_reversed_indices
from src.exon import Exon

class TestExon(unittest.TestCase):

    def setUp(self):
        test_identifier1 = 3
        test_indices1 = [3734, 4034] 
        test_score1 = 0.9
        test_parent_id1 = 2
        self.test_exon0 = Exon(identifier=test_identifier1, indices=test_indices1, score=test_score1, parent_id=test_parent_id1)
        self.extra_identifiers = [4, 5, 6, 7]
        self.extra_scores = [0.9, 0.9, 0.9, 0.9]
        self.extra_indices = [[4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        self.test_exon1 = Exon(identifier=test_identifier1, indices=test_indices1, score=test_score1, parent_id=test_parent_id1)
        for ind_pair in self.extra_indices:
            self.test_exon1.add_indices(ind_pair)
        for ident in self.extra_identifiers:
            self.test_exon1.add_identifier(ident)
        for score in self.extra_scores:
            self.test_exon1.add_score(score)

    def test_constructor(self):
        self.assertEquals('Exon', self.test_exon1.__class__.__name__)
        self.assertEquals('+', self.test_exon1.strand)

    def test_add_indices(self):
        for ind_pair in self.extra_indices:
            self.test_exon0.add_indices(ind_pair)
        self.assertEquals(5, len(self.test_exon0.indices))
        self.assertEquals([5249, 6565], self.test_exon0.indices[3])

    def test_add_identifier(self):
        for ident in self.extra_identifiers:
            self.test_exon0.add_identifier(ident)
        self.assertEquals(5, len(self.test_exon0.identifier))
        self.assertEquals(6, self.test_exon0.identifier[3])

    def test_add_score(self):
        for score in self.extra_scores:
            self.test_exon0.add_score(score)
        self.assertEquals(5, len(self.test_exon0.score))
        self.assertEquals(0.9, self.test_exon0.score[4])
       
    def test_length(self): 
        self.assertEquals(301, self.test_exon0.length())
        self.assertEquals(3453, self.test_exon1.length())
       
    def test_adjust_indices(self): 
        self.test_exon1.adjust_indices(7)
        self.assertEquals(4339, self.test_exon1.indices[1][1])
        # adjust them back...
        self.test_exon1.adjust_indices(-7)
        self.assertEquals(4332, self.test_exon1.indices[1][1])

    def test_to_gff(self):
        expected1 = "sctg_0080_0020\tmaker\texon\t3734\t4034\t0.9\t+\t.\tID=3;Parent=2\n"
        expected2 = "sctg_0080_0020\tmaker\texon\t4092\t4332\t0.9\t+\t.\tID=4;Parent=2\n"
        expected3 = "sctg_0080_0020\tmaker\texon\t4399\t5185\t0.9\t+\t.\tID=5;Parent=2\n"
        expected4 = "sctg_0080_0020\tmaker\texon\t5249\t6565\t0.9\t+\t.\tID=6;Parent=2\n"
        expected5 = "sctg_0080_0020\tmaker\texon\t6630\t7436\t0.9\t+\t.\tID=7;Parent=2\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = self.test_exon1.to_gff(seq_name="sctg_0080_0020", source="maker")
        self.assertEquals(expected, actual)

    def test_to_tbl_positive_complete(self):
        expected = "3734\t4034\tmRNA\n"
        expected += "4092\t4332\n"
        expected += "4399\t5185\n"
        expected += "5249\t6565\n"
        expected += "6630\t7436\n"
        self.test_exon1.strand = '+'
        self.assertEquals(self.test_exon1.to_tbl(True, True, "mRNA"), expected)

    def test_to_tbl_positive_no_start(self):
        expected = "<3734\t4034\tmRNA\n"
        expected += "4092\t4332\n"
        expected += "4399\t5185\n"
        expected += "5249\t6565\n"
        expected += "6630\t7436\n"
        self.test_exon1.strand = '+'
        self.assertEquals(self.test_exon1.to_tbl(False, True, "mRNA"), expected)

    def test_to_tbl_negative_complete(self):
        expected = "7436\t6630\tmRNA\n"
        expected += "6565\t5249\n"
        expected += "5185\t4399\n"
        expected += "4332\t4092\n"
        expected += "4034\t3734\n"
        self.test_exon1.strand = '-'
        self.assertEquals(self.test_exon1.to_tbl(True, True, "mRNA"), expected)

    def test_to_tbl_negative_no_start_no_stop(self):
        expected = "<7436\t6630\tmRNA\n"
        expected += "6565\t5249\n"
        expected += "5185\t4399\n"
        expected += "4332\t4092\n"
        expected += "4034\t>3734\n"
        self.test_exon1.strand = '-'
        self.assertEquals(self.test_exon1.to_tbl(False, False, "mRNA"), expected)

    def test_get_reversed_indices(self):
        indices = [[1, 10], [20, 30], [40, 50]]
        expected = [[50, 40], [30, 20], [10, 1]]
        self.assertEquals(get_reversed_indices(indices), expected)

        

##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestExon))
    return suite

if __name__ == '__main__':
    unittest.main()
