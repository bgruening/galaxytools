#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.cds import CDS

class TestCDS(unittest.TestCase):

    def setUp(self):
        self.test_indices1 = [3734, 4034]
        self.extra_indices = [[4092, 4332], [4399, 5185], [5249, 6565], [6630, 7436]]
        test_identifier1 = 8
        self.extra_identifiers = [9, 10, 11, 12]
        self.test_phase1 = 0
        self.extra_phases = [2, 1, 0, 0]
        test_parent_id1 = 2
        self.test_cds0 = CDS(identifier=test_identifier1, indices=self.test_indices1, score=None, phase=self.test_phase1, strand='-', parent_id=test_parent_id1)
        self.test_cds1 = CDS(identifier=test_identifier1, indices=self.test_indices1, score=None, phase=self.test_phase1, strand='+', parent_id=test_parent_id1)
        for ind_pair in self.extra_indices:
            self.test_cds1.add_indices(ind_pair)
        for ident in self.extra_identifiers:
            self.test_cds1.add_identifier(ident)
        for phase in self.extra_phases:
            self.test_cds1.add_phase(phase)

    def test_constructor(self):
        self.assertEquals('-', self.test_cds0.strand)

    def test_get_start_indices_pos_strand(self):
        expected = [3734, 3736]
        self.assertEquals(expected, self.test_cds1.get_start_indices('+'))

    def test_get_start_indices_neg_strand(self):
        expected = [4032, 4034]
        self.assertEquals(expected, self.test_cds1.get_start_indices('-'))

    def test_get_stop_indices_pos_strand(self):
        expected = [7434, 7436]
        self.assertEquals(expected, self.test_cds1.get_stop_indices('+'))

    def test_get_stop_indices_neg_strand(self):
        expected = [6630, 6632]
        self.assertEquals(expected, self.test_cds1.get_stop_indices('-'))

    def test_extract_sequence_pos_strand(self):
        seq_object = Mock()
        seq_object.get_subseq.return_value = 'GATTACA'
        strand = '+'
        seq = self.test_cds1.extract_sequence(seq_object, strand)
        expected = 'GATTACAGATTACAGATTACAGATTACAGATTACA'
        self.assertEquals(expected, seq)

    def test_extract_sequence_neg_strand(self):
        seq_object = Mock()
        seq_object.get_subseq.return_value = 'GATTACA'
        strand = '-'
        result = self.test_cds1.extract_sequence(seq_object, strand)
        expected = 'TGTAATCTGTAATCTGTAATCTGTAATCTGTAATC'
        self.assertEquals(expected, result)

    def test_cds_constructor(self):
        self.assertEquals('CDS', self.test_cds0.__class__.__name__)
        # should also be able to construct w/o all the params...
        empty_cds = CDS()
        self.assertEquals('CDS', empty_cds.feature_type) 

    def test_add_indices(self):
        for ind_pair in self.extra_indices:
            self.test_cds0.add_indices(ind_pair)
        self.assertEquals([4399, 5185], self.test_cds0.indices[2])

    def test_add_identifier(self):
        for ident in self.extra_identifiers:
            self.test_cds0.add_identifier(ident)
        self.assertEquals(5, len(self.test_cds0.identifier))

    def test_add_phase(self):
        for phase in self.extra_phases:
            self.test_cds0.add_phase(phase)
        self.assertEquals(5, len(self.test_cds0.phase))
        self.assertEquals(1, self.test_cds0.phase[2])

    def test_sort_attributes(self):
        cds = CDS()
        cds.indices = [[25, 30], [5, 10]] # out of order!
        cds.identifier = ["cds2", "cds1"]
        cds.phase = [1, 0]
        self.assertEquals("cds1", cds.identifier[1])
        self.assertEquals([25, 30], cds.indices[0])
        self.assertEquals(1, cds.phase[0])
        cds.sort_attributes()
        self.assertEquals("cds1", cds.identifier[0])
        self.assertEquals([5, 10], cds.indices[0])
        self.assertEquals(0, cds.phase[0])

    def test_length(self):
        self.assertEquals(3453, self.test_cds1.length())

    def test_adjust_indices(self):
        self.test_cds1.adjust_indices(-5)
        self.assertEquals(3729, self.test_cds1.indices[0][0])
        # (adjust them back so future test don't get confused :)
        self.test_cds1.adjust_indices(5)
        self.assertEquals(5185, self.test_cds1.indices[2][1])

    def test_to_gff(self):
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=8;Parent=2;foo=dog\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=9;Parent=2;foo=dog\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=10;Parent=2;foo=dog\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=11;Parent=2;foo=dog\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=12;Parent=2;foo=dog\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        self.test_cds1.add_annotation('foo','dog') # Make sure our annotations are working
        actual = self.test_cds1.to_gff(seq_name="sctg_0080_0020", source="maker")
        self.assertEquals(expected, actual)
        # what if identifier, parent_id are strings? does it matter?
        test_cds2 = CDS(identifier='foo1', indices=self.test_indices1, score=None, strand='+', phase=self.test_phase1, parent_id='bar7')
        extra_identifiers2 = ['foo2', 'foo3', 'foo4', 'foo5']
        for ind_pair in self.extra_indices:
            test_cds2.add_indices(ind_pair)
        for ident in extra_identifiers2:
            test_cds2.add_identifier(ident)
        for phase in self.extra_phases:
            test_cds2.add_phase(phase)
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Parent=bar7\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Parent=bar7\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Parent=bar7\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Parent=bar7\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Parent=bar7\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker")
        self.assertEquals(expected, actual)
        expected1 = "sctg_0080_0020\tmaker\tCDS\t3734\t4034\t.\t+\t0\tID=foo1;Parent=bar7\n"
        expected2 = "sctg_0080_0020\tmaker\tCDS\t4092\t4332\t.\t+\t2\tID=foo2;Parent=bar7\n"
        expected3 = "sctg_0080_0020\tmaker\tCDS\t4399\t5185\t.\t+\t1\tID=foo3;Parent=bar7\n"
        expected4 = "sctg_0080_0020\tmaker\tCDS\t5249\t6565\t.\t+\t0\tID=foo4;Parent=bar7\n"
        expected5 = "sctg_0080_0020\tmaker\tCDS\t6630\t7436\t.\t+\t0\tID=foo5;Parent=bar7\n"
        expected = expected1 + expected2 + expected3 + expected4 + expected5
        actual = test_cds2.to_gff(seq_name="sctg_0080_0020", source="maker")
        self.assertEquals(expected, actual)

    def test_to_tbl_positive_complete(self):
        expected = "3734\t4034\tCDS\n"
        expected += "4092\t4332\n"
        expected += "4399\t5185\n"
        expected += "5249\t6565\n"
        expected += "6630\t7436\n"
        expected += "\t\t\tcodon_start\t2\n"
        self.test_cds1.phase[0] = 1
        self.assertEquals(self.test_cds1.to_tbl(True, True), expected)

    def test_to_tbl_negative_complete(self):
        expected = "7436\t6630\tCDS\n"
        expected += "6565\t5249\n"
        expected += "5185\t4399\n"
        expected += "4332\t4092\n"
        expected += "4034\t3734\n"
        expected += "\t\t\tcodon_start\t1\n"
        self.test_cds1.strand = '-'
        self.assertEquals(self.test_cds1.to_tbl(True, True), expected)

    def test_to_tbl_negative_no_start_no_stop(self):
        expected = "<7436\t6630\tCDS\n"
        expected += "6565\t5249\n"
        expected += "5185\t4399\n"
        expected += "4332\t4092\n"
        expected += "4034\t>3734\n"
        expected += "\t\t\tcodon_start\t2\n"
        # shouldn't look at phase[0] for negative strand!
        self.test_cds1.phase[0] = 2 # should ignore this.
        self.test_cds1.phase[4] = 1
        self.test_cds1.strand = '-'
        self.assertEquals(self.test_cds1.to_tbl(False, False), expected)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestCDS))
    return suite

if __name__ == '__main__':
    unittest.main()
