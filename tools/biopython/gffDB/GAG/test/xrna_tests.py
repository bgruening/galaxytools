#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.xrna import XRNA

class TestXRNA(unittest.TestCase):

    def setUp(self):
        self.test_mrna0 = XRNA(identifier='bdor_foo', indices=[3734, 7436], strand='-', parent_id=1)
        self.test_mrna1 = XRNA(identifier='bdor_foo2', indices=[3734, 7436], parent_id=1, seq_name="sctg_0080_0020", source="maker")
        self.fake_exon = Mock()
        self.fake_cds = Mock()
        self.fake_start_codon = Mock()
        self.test_mrna1.exon = self.fake_exon
        self.test_mrna1.cds = self.fake_cds
        self.test_mrna1.add_other_feature(self.fake_start_codon)

    def test_constructor(self):
        self.assertEquals('XRNA', self.test_mrna0.__class__.__name__)
        self.assertEquals('-', self.test_mrna0.strand)

    def test_add_annotation(self):
        expected = {'foo': ['dog', 'cat']}
        self.test_mrna0.add_annotation('foo', 'dog')
        self.test_mrna0.add_annotation('foo', 'cat')
        self.assertEqual(expected, self.test_mrna0.annotations)

    def test_length(self):
        self.assertEqual(3703, self.test_mrna0.length())

    def test_number_of_gagflags(self):
        self.fake_cds.gagflagged.return_value = False
        self.fake_exon.gagflagged.return_value = True
        self.assertEquals(1, self.test_mrna1.number_of_gagflags())

    def test_add_other_feature(self): 
        self.assertEquals(0, len(self.test_mrna0.other_features))
        self.test_mrna0.add_other_feature(self.fake_start_codon)
        self.assertEquals(1, len(self.test_mrna0.other_features))

    def test_add_start_codon(self):
        self.assertFalse(self.test_mrna0.has_start())
        self.test_mrna0.add_start_codon([4000, 4002])
        self.assertTrue(self.test_mrna0.has_start())
        self.assertEquals([[4000, 4002]], self.test_mrna0.other_features[0].indices)

    def test_add_start_codon_gets_the_strand_right(self):
        self.assertFalse(self.test_mrna0.has_start())
        self.test_mrna0.add_start_codon([4000, 4002])
        self.assertTrue(self.test_mrna0.has_start())
        self.assertEquals([[4000, 4002]], self.test_mrna0.other_features[0].indices)
        self.assertEquals('-', self.test_mrna0.other_features[0].strand)

    def test_add_stop_codon(self):
        self.assertFalse(self.test_mrna0.has_stop())
        self.test_mrna0.add_stop_codon([7000, 7002])
        self.assertTrue(self.test_mrna0.has_stop())

    def test_adjust_indices(self):
        self.test_mrna1.adjust_indices(32)
        self.assertEqual(7468, self.test_mrna1.indices[1])
        self.fake_cds.adjust_indices.assert_called_with(32, 1)
        self.fake_exon.adjust_indices.assert_called_with(32, 1)
        self.fake_start_codon.adjust_indices.assert_called_with(32, 1)
        # adjust 'em back
        self.test_mrna1.adjust_indices(-32)
        self.assertEqual(3734, self.test_mrna1.indices[0])
        self.fake_cds.adjust_indices.assert_called_with(-32, 1)
        self.fake_exon.adjust_indices.assert_called_with(-32, 1)
        self.fake_start_codon.adjust_indices.assert_called_with(-32, 1)

    def test_adjust_indices_after_start_index(self):
        self.test_mrna1.adjust_indices(32, 8000)
        self.assertEqual(3734, self.test_mrna1.indices[0])
        self.test_mrna1.adjust_indices(32, 3000)
        self.assertEqual(3766, self.test_mrna1.indices[0])

    def test_has_start(self):
        self.fake_start_codon.feature_type = 'start_codon'
        self.assertTrue(self.test_mrna1.has_start())

    def test_has_stop(self):
        self.assertFalse(self.test_mrna1.has_stop())

    def test_str(self):
        expected = "mRNA (ID=bdor_foo2) containing Exon, CDS and 1 other features"
        self.assertEquals(expected, str(self.test_mrna1))

    def test_cds_to_gff(self):
        foo = self.test_mrna1.cds_to_gff("foo_seq", "maker")
        self.fake_cds.to_gff.assert_called_with("foo_seq", "maker")

    def test_cds_to_tbl(self):
        foo = self.test_mrna1.cds_to_tbl()
        self.fake_cds.to_tbl.assert_called_with(False, False)

    def test_to_gff(self):
        self.fake_exon.to_gff.return_value = "...exon to gff\n"
        self.fake_cds.to_gff.return_value = "...cds to gff\n"
        self.fake_start_codon.to_gff.return_value = "...start codon to gff\n"
        expected = "sctg_0080_0020\tmaker\tmRNA\t"
        expected += "3734\t7436\t.\t+\t.\t"
        expected += "ID=bdor_foo2;Parent=1;foo=dog\n"
        expected += "...exon to gff\n...cds to gff\n"
        expected += "...start codon to gff\n"
        self.test_mrna1.add_annotation('foo', 'dog')
        actual = self.test_mrna1.to_gff()
        self.assertEquals(expected, actual)
        self.fake_exon.to_gff.assert_called_with("sctg_0080_0020", "maker")
        self.fake_cds.to_gff.assert_called_with("sctg_0080_0020", "maker")
        self.fake_start_codon.to_gff.assert_called_with("sctg_0080_0020", "maker")

    def test_to_gff_multiple_dbxref(self):
        self.fake_exon.to_gff.return_value = "...exon to gff\n"
        self.fake_cds.to_gff.return_value = "...cds to gff\n"
        self.fake_start_codon.to_gff.return_value = "...start codon to gff\n"
        expected = "sctg_0080_0020\tmaker\tmRNA\t"
        expected += "3734\t7436\t.\t+\t.\t"
        expected += "ID=bdor_foo2;Parent=1;foo=dog,cat\n"
        expected += "...exon to gff\n...cds to gff\n"
        expected += "...start codon to gff\n"
        self.test_mrna1.add_annotation('foo', 'dog')
        self.test_mrna1.add_annotation('foo', 'cat')
        actual = self.test_mrna1.to_gff()
        self.assertEquals(expected, actual)
        self.fake_exon.to_gff.assert_called_with("sctg_0080_0020", "maker")
        self.fake_cds.to_gff.assert_called_with("sctg_0080_0020", "maker")
        self.fake_start_codon.to_gff.assert_called_with("sctg_0080_0020", "maker")


    def test_indices_intersect_mrna_false(self):
        mrna = XRNA(identifier=1, indices=[10, 20], parent_id='foo')
        self.assertFalse(mrna.indices_intersect_mrna([5, 9]))
        self.assertFalse(mrna.indices_intersect_mrna([21, 25]))

    def test_indices_intersect_mrna_true(self):
        mrna = XRNA(identifier=1, indices=[10, 20], parent_id='foo')
        self.assertTrue(mrna.indices_intersect_mrna([5, 10]))
        self.assertTrue(mrna.indices_intersect_mrna([20, 25]))
        self.assertTrue(mrna.indices_intersect_mrna([9, 21]))

    def test_create_start_and_stop_if_necessary(self):
        seq_object = Mock()
        cds = Mock()
        cds.extract_sequence.return_value = 'atgtag' # startstop
        cds.get_start_indices.return_value = 20
        cds.get_stop_indices.return_value = 40
        self.test_mrna0.cds = cds
        seq_name = 'seq_foo'
        strand = '+'
        self.assertFalse(self.test_mrna0.other_features)
        self.test_mrna0.create_start_and_stop_if_necessary(seq_object, strand)
        self.assertTrue(self.test_mrna0.other_features)
        self.assertEquals(2, len(self.test_mrna0.other_features))

    def test_create_start_and_stop_when_no_start_or_stop(self):
        seq_object = Mock()
        cds = Mock()
        cds.extract_sequence.return_value = 'tagatg' # no start or stop
        cds.get_start_indices.return_value = 20
        cds.get_stop_indices.return_value = 40
        self.test_mrna0.cds = cds
        seq_name = 'seq_foo'
        strand = '+'
        self.assertFalse(self.test_mrna0.other_features)
        self.test_mrna0.create_start_and_stop_if_necessary(seq_object, strand)
        self.assertFalse(self.test_mrna0.other_features)

    def test_to_tbl(self):
        self.fake_exon.to_tbl.return_value = "fake_exon_to_tbl...\n"
        self.fake_cds.to_tbl.return_value = "fake_cds_to_tbl...\n"
        expected = "fake_exon_to_tbl...\n"
        expected += "\t\t\tproduct\thypothetical protein\n"
        expected += "\t\t\tprotein_id\tgnl|ncbi|bdor_foo2\n"
        expected += "\t\t\ttranscript_id\tgnl|ncbi|bdor_foo2_mrna\n"
        expected += "fake_cds_to_tbl...\n"
        expected += "\t\t\tproduct\thypothetical protein\n"
        expected += "\t\t\tprotein_id\tgnl|ncbi|bdor_foo2\n"
        expected += "\t\t\ttranscript_id\tgnl|ncbi|bdor_foo2_mrna\n"
        self.assertEquals(self.test_mrna1.to_tbl(), expected)

    def test_to_tbl_with_annotations(self):
        self.fake_exon.to_tbl.return_value = "fake_exon_to_tbl...\n"
        self.fake_cds.to_tbl.return_value = "fake_cds_to_tbl...\n"
        self.test_mrna1.add_annotation('foo', 'dog')
        self.test_mrna1.add_annotation('foo', 'cat')
        expected = "fake_exon_to_tbl...\n"
        expected += "\t\t\tproduct\thypothetical protein\n"
        expected += "\t\t\tprotein_id\tgnl|ncbi|bdor_foo2\n"
        expected += "\t\t\ttranscript_id\tgnl|ncbi|bdor_foo2_mrna\n"
        expected += "fake_cds_to_tbl...\n"
        expected += "\t\t\tfoo\tdog\n"
        expected += "\t\t\tfoo\tcat\n"
        expected += "\t\t\tproduct\thypothetical protein\n"
        expected += "\t\t\tprotein_id\tgnl|ncbi|bdor_foo2\n"
        expected += "\t\t\ttranscript_id\tgnl|ncbi|bdor_foo2_mrna\n"
        self.assertEquals(self.test_mrna1.to_tbl(), expected)

    def test_to_tbl_with_product(self):
        self.fake_exon.to_tbl.return_value = "fake_exon_to_tbl...\n"
        self.fake_cds.to_tbl.return_value = "fake_cds_to_tbl...\n"
        self.test_mrna1.add_annotation('product', 'dog')
        expected = "fake_exon_to_tbl...\n"
        expected += "\t\t\tproduct\tdog\n"
        expected += "\t\t\tprotein_id\tgnl|ncbi|bdor_foo2\n"
        expected += "\t\t\ttranscript_id\tgnl|ncbi|bdor_foo2_mrna\n"
        expected += "fake_cds_to_tbl...\n"
        expected += "\t\t\tproduct\tdog\n"
        expected += "\t\t\tprotein_id\tgnl|ncbi|bdor_foo2\n"
        expected += "\t\t\ttranscript_id\tgnl|ncbi|bdor_foo2_mrna\n"
        print(expected)
        print(self.test_mrna1.to_tbl())
        self.assertEquals(expected, self.test_mrna1.to_tbl())

    def test_indices_intersect_cds_false(self):
        self.fake_cds.indices_intersect_cds.return_value = False
        self.assertFalse(self.test_mrna1.indices_intersect_cds([1, 9]))

    def test_indices_intersect_cds_true(self):
        self.fake_cds.indices_intersect_cds.return_value = True
        self.assertTrue(self.test_mrna1.indices_intersect_cds([1, 9]))

    def test_annotations_contains_product(self):
        self.assertFalse(self.test_mrna0.annotations_contain_product())
        self.test_mrna0.add_annotation('product', 'foo')
        self.assertTrue(self.test_mrna0.annotations_contain_product())


    ## STATS STUFF ##

    def set_fake_exon_indices(self):
        self.fake_exon.indices = [[1, 5], [11, 16], [21, 27]]

    def test_get_longest_exon(self):
        self.set_fake_exon_indices()
        self.assertEquals(7, self.test_mrna1.get_longest_exon())

    def test_get_shortest_exon(self):
        self.set_fake_exon_indices()
        self.assertEquals(5, self.test_mrna1.get_shortest_exon())

    def test_get_total_exon_length(self):
        self.set_fake_exon_indices()
        self.assertEquals(18, self.test_mrna1.get_total_exon_length())

    def test_get_num_exons(self):
        self.set_fake_exon_indices()
        self.assertEquals(3, self.test_mrna1.get_num_exons())

    def test_get_num_exons_if_no_exons(self):
        # test_mrna0 has no exon
        self.assertEquals(0, self.test_mrna0.get_num_exons())

    def test_nothing_fails_if_no_exons(self):
        self.assertFalse(self.test_mrna0.get_longest_exon())
        self.assertFalse(self.test_mrna0.get_shortest_exon())
        self.assertFalse(self.test_mrna0.get_total_exon_length())
        self.assertFalse(self.test_mrna0.get_longest_intron())
        self.assertFalse(self.test_mrna0.get_shortest_intron())
        self.assertFalse(self.test_mrna0.get_total_intron_length())

    def test_get_longest_intron(self):
        self.set_fake_exon_indices()
        self.assertEquals(5, self.test_mrna1.get_longest_intron())

    def test_get_shortest_intron(self):
        self.set_fake_exon_indices()
        self.assertEquals(4, self.test_mrna1.get_shortest_intron())

    def test_get_total_intron_length(self):
        self.set_fake_exon_indices()
        self.assertEquals(13, self.test_mrna1.get_total_intron_length())

    def test_get_num_introns(self):
        self.set_fake_exon_indices()
        self.assertEquals(2, self.test_mrna1.get_num_introns())



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestXRNA))
    return suite

if __name__ == '__main__':
    unittest.main()
