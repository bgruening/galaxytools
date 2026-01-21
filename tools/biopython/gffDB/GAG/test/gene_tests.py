#!/usr/bin/env python

import unittest
from mock import Mock, PropertyMock
from src.gene import Gene

class TestGene(unittest.TestCase):

    def setUp(self):
        self.test_gene0 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1)
        self.test_gene1 = Gene(seq_name="sctg_0080_0020", source="maker", indices=[3734, 7436], strand='+', identifier=1)
        
        self.fake_mrna1 = Mock()
        self.fake_mrna1.identifier = "fake_mrna1"
        self.fake_mrna1.death_flagged = False
        
        self.fake_mrna2 = Mock()
        self.fake_mrna2.identifier = "fake_mrna2"
        self.fake_mrna2.death_flagged = False
        
        self.test_gene1.mrnas.append(self.fake_mrna1)
        self.test_gene1.mrnas.append(self.fake_mrna2)

    def test_constructor(self):
        self.assertEqual('Gene', self.test_gene0.__class__.__name__)

    def test_length(self):
        self.assertEqual(3703, self.test_gene0.length())

    def test_gagflagged(self):
        self.assertFalse(self.test_gene0.gagflagged())
        self.test_gene0.annotations = {"gag_flag": ["nice gene"]}
        self.assertTrue(self.test_gene0.gagflagged())

    def test_number_of_gagflags(self):
        self.fake_mrna1.number_of_gagflags.return_value = 2
        self.fake_mrna2.number_of_gagflags.return_value = 1
        self.test_gene1.annotations = {"gag_flag": ["nice gene"]}
        self.assertEquals(4, self.test_gene1.number_of_gagflags())

    def test_get_mrna_ids(self):
        expected = ["fake_mrna1", "fake_mrna2"]
        self.assertEquals(self.test_gene1.get_mrna_ids(), expected)
    
    def test_remove_mrna(self):
        self.assertEquals(self.test_gene1.mrnas, [self.fake_mrna1, self.fake_mrna2])
        self.assertEquals(len(self.test_gene1.removed_mrnas), 0)
        self.test_gene1.remove_mrna('fake_mrna1')
        self.assertEquals(self.test_gene1.mrnas, [self.fake_mrna2])
        self.assertEquals(self.test_gene1.removed_mrnas, [self.fake_mrna1])

    def test_remove_mrnas_from_list(self):
        self.fake_mrna3 = Mock()
        self.fake_mrna3.identifier = "fake_mrna3"
        self.fake_mrna3.death_flagged = False
        self.test_gene1.mrnas.append(self.fake_mrna3)
        bad_mrnas = ["fake_mrna3", "fake_mrna1"]
        self.assertEquals(3, len(self.test_gene1.mrnas))
        self.assertEquals(0, len(self.test_gene1.removed_mrnas))
        removed_mrnas = self.test_gene1.remove_mrnas_from_list(bad_mrnas)
        self.assertEquals(2, len(removed_mrnas))
        self.assertEquals(1, len(self.test_gene1.mrnas))
        self.assertEquals(2, len(self.test_gene1.removed_mrnas))
    
    def test_remove_empty_mrnas(self):
        self.fake_mrna1.rna_type = "mRNA"
        self.fake_mrna1.cds = Mock()
        self.fake_mrna1.exon = None
        self.fake_mrna2.rna_type = "mRNA"
        self.fake_mrna2.cds = None
        self.fake_mrna2.exon = Mock()
        self.fake_mrna3 = Mock()
        self.fake_mrna3.rna_type = "mRNA"
        self.fake_mrna3.identifier = "fake_mrna3"
        self.fake_mrna3.death_flagged = False
        self.fake_mrna3.cds = Mock()
        self.fake_mrna3.exon = Mock()
        self.test_gene1.mrnas.append(self.fake_mrna3)
        self.assertEquals(3, len(self.test_gene1.mrnas))
        self.assertEquals(0, len(self.test_gene1.removed_mrnas))
        removed_mrnas = self.test_gene1.remove_empty_mrnas()
        self.assertEquals(2, len(removed_mrnas))
        self.assertEquals(1, len(self.test_gene1.mrnas))
        self.assertEquals(2, len(self.test_gene1.removed_mrnas))

    def test_get_longest_exon(self):
        self.fake_mrna1.get_longest_exon.return_value = 10
        self.fake_mrna2.get_longest_exon.return_value = 20
        self.assertEquals(20, self.test_gene1.get_longest_exon())

    def test_get_shortest_exon(self):
        self.fake_mrna1.get_shortest_exon.return_value = 5
        self.fake_mrna2.get_shortest_exon.return_value = 8
        self.assertEquals(5, self.test_gene1.get_shortest_exon())

    def test_get_total_exon_length(self):
        self.fake_mrna1.get_total_exon_length.return_value = 15
        self.fake_mrna2.get_total_exon_length.return_value = 25
        self.assertEquals(40, self.test_gene1.get_total_exon_length())

    def test_get_num_exons(self):
        self.fake_mrna1.get_num_exons.return_value = 5
        self.fake_mrna2.get_num_exons.return_value = 4
        self.assertEquals(9, self.test_gene1.get_num_exons())

    def test_get_longest_intron(self):
        self.fake_mrna1.get_longest_intron.return_value = 8
        self.fake_mrna2.get_longest_intron.return_value = 10
        self.assertEquals(10, self.test_gene1.get_longest_intron())

    def test_get_shortest_intron(self):
        self.fake_mrna1.get_shortest_intron.return_value = 5
        self.fake_mrna2.get_shortest_intron.return_value = 8
        self.assertEquals(5, self.test_gene1.get_shortest_intron())

    def test_get_total_intron_length(self):
        self.fake_mrna1.get_total_intron_length.return_value = 15
        self.fake_mrna2.get_total_intron_length.return_value = 25
        self.assertEquals(40, self.test_gene1.get_total_intron_length())

    def test_get_num_introns(self):
        self.fake_mrna1.get_num_introns.return_value = 3
        self.fake_mrna2.get_num_introns.return_value = 2
        self.assertEquals(5, self.test_gene1.get_num_introns())

    def test_get_partial_info(self):
        self.fake_mrna1.has_stop.return_value = True
        self.fake_mrna1.has_start.return_value = True
        self.fake_mrna2.has_stop.return_value = False
        self.fake_mrna2.has_start.return_value = True
        results = self.test_gene1.get_partial_info()
        self.assertEquals(1, results["complete"])

    def test_adjust_indices(self):
        self.test_gene1.adjust_indices(16)
        self.fake_mrna1.adjust_indices.assert_called_with(16, 1)
        self.assertEquals(3750, self.test_gene1.indices[0])
        # adjust them back
        self.test_gene1.adjust_indices(-16)
        self.fake_mrna1.adjust_indices.assert_called_with(-16, 1)
        self.assertEquals(3734, self.test_gene1.indices[0])

    def test_remove_mrnas_with_internal_stops(self):
        helper = Mock()
        helper.mrna_contains_internal_stop.return_value = True
        self.assertEquals(2, len(self.test_gene1.mrnas))
        self.test_gene1.remove_mrnas_with_internal_stops(helper)
        self.assertEquals(0, len(self.test_gene1.mrnas))

    def test_contains_mrna(self):
        self.fake_mrna1.identifier = "foo_mrna"
        self.fake_mrna2.identifier = "bar_mrna"
        self.assertTrue(self.test_gene1.contains_mrna("foo_mrna"))
        self.assertFalse(self.test_gene1.contains_mrna("zub_mrna"))

    def test_cds_to_gff(self):
        self.fake_mrna1.identifier = "foo_mrna"
        foo = self.test_gene1.cds_to_gff("foo_seq", "foo_mrna")
        self.fake_mrna1.cds_to_gff.assert_called_with("foo_seq", "maker")

    def test_cds_to_gff_no_such_mrna(self):
        self.fake_mrna1.identifier = "foo_mrna"
        foo = self.test_gene1.cds_to_gff("foo_seq", "bar_mrna")
        self.assertFalse(foo)

    def test_cds_to_tbl(self):
        self.fake_mrna1.identifier = "foo_mrna"
        foo = self.test_gene1.cds_to_tbl("foo_mrna")
        self.fake_mrna1.cds_to_tbl.assert_called_with()

    def test_to_mrna_fasta(self):
        helper = Mock()
        helper.mrna_to_fasta.return_value = "mrna_to_fasta\n"
        expected = "mrna_to_fasta\nmrna_to_fasta\n"
        self.assertEquals(expected, self.test_gene1.to_mrna_fasta(helper))

    def test_to_cds_fasta(self):
        helper = Mock()
        helper.mrna_to_cds_fasta.return_value = "mrna_to_CDS_fasta\n"
        expected = "mrna_to_CDS_fasta\nmrna_to_CDS_fasta\n"
        self.assertEquals(expected, self.test_gene1.to_cds_fasta(helper))

    def test_to_protein_fasta(self):
        helper = Mock()
        helper.mrna_to_protein_fasta.return_value = "mrna_to_protein_fasta\n"
        expected = "mrna_to_protein_fasta\nmrna_to_protein_fasta\n"
        self.assertEquals(expected, self.test_gene1.to_protein_fasta(helper))

    def test_to_gff(self):
        self.fake_mrna1.to_gff.return_value = "fake mrna1 to gff here:)\n"
        self.fake_mrna2.to_gff.return_value = "fake mrna2 to gff here:)\n"
        expected = "sctg_0080_0020\tmaker\tgene\t3734\t7436\t.\t+\t."
        expected += "\tID=1;foo=dog\n"
        expected += "fake mrna1 to gff here:)\n"
        expected += "fake mrna2 to gff here:)\n"
        self.test_gene1.add_annotation('foo', 'dog')
        self.assertEquals(expected, self.test_gene1.to_gff())

    def test_to_gff_with_name(self):
        self.fake_mrna1.to_gff.return_value = "fake mrna1 to gff here:)\n"
        self.fake_mrna2.to_gff.return_value = "fake mrna2 to gff here:)\n"
        expected = "sctg_0080_0020\tmaker\tgene\t3734\t7436\t.\t+\t."
        expected += "\tID=1;Name=foo_gene;foo=dog\n"
        expected += "fake mrna1 to gff here:)\n"
        expected += "fake mrna2 to gff here:)\n"
        self.test_gene1.add_annotation('foo', 'dog')
        self.test_gene1.name = "foo_gene"
        self.assertEquals(expected, self.test_gene1.to_gff())

    def test_str(self):
        expected = "Gene (ID=1, seq_name=sctg_0080_0020) containing 2 mrnas"
        self.assertEquals(expected, str(self.test_gene1))

    def test_create_starts_and_stops(self):
        mrna1 = Mock()
        mrna2 = Mock()
        self.test_gene0.mrnas = [mrna1, mrna2]
        seq_object = Mock()
        self.test_gene0.create_starts_and_stops(seq_object)
        mrna1.create_start_and_stop_if_necessary.assert_called_with(seq_object, '+')
        mrna2.create_start_and_stop_if_necessary.assert_called_with(seq_object, '+')

    def test_add_mrna_annotation(self):
        mrna = Mock()
        mrna.identifier = "foo_mrna"
        self.test_gene0.mrnas = [mrna]
        self.test_gene0.add_mrna_annotation("foo_mrna", "gag_flag", "awesome_anno")
        mrna.add_annotation.assert_called_with("gag_flag", "awesome_anno")

    def test_to_tbl_positive(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "1\t50\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

    def test_to_tbl_positive_start_nostop(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        mrna2.has_start.return_value = True
        mrna2.has_stop.return_value = False
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "1\t>50\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

    def test_to_tbl_positive_nostart_stop(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        mrna2.has_start.return_value = False
        mrna2.has_stop.return_value = True
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "<1\t50\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

    def test_to_tbl_positive_nostart_nostop(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        mrna2.has_start.return_value = False
        mrna2.has_stop.return_value = False
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "<1\t>50\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)


    def test_to_tbl_positive_with_name(self):
        gene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1", name="wtfg")
        self.assertFalse(gene.annotations)
        gene.add_annotation('foo', 'dog')
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "1\t50\tgene\n\t\t\tgene\twtfg\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

    def test_to_tbl_negative(self):
        gene = Gene("seq1", "maker", [1, 50], "-", "foo_gene_1")
        mrna1 = Mock()
        mrna1.to_tbl.return_value = "mrna1_to_tbl...\n"
        mrna2 = Mock()
        mrna2.to_tbl.return_value = "mrna2_to_tbl...\n"
        gene.mrnas.append(mrna1)
        gene.mrnas.append(mrna2)
        expected = "50\t1\tgene\n\t\t\tlocus_tag\tfoo_gene_1\nmrna1_to_tbl...\nmrna2_to_tbl...\n"
        self.assertEquals(gene.to_tbl(), expected)

    def test_gene_initialized_without_annotations(self):
        newgene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1")
        self.assertFalse(newgene.annotations)
        self.assertEquals(0, len(newgene.annotations.keys()))

    def test_gene_initialized_with_annotations(self):
        newgene = Gene(seq_name="seq1", source="maker", indices=[1, 50], strand="+", identifier="foo_gene_1",\
                annotations={"bar": ["cat"]})
        self.assertTrue(newgene.annotations)
        self.assertEquals(1, len(newgene.annotations.keys()))


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGene))
    return suite

if __name__ == '__main__':
    unittest.main()
