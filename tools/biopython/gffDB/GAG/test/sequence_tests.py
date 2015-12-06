#!/usr/bin/env python

import unittest
from mock import Mock, MagicMock
from src.sequence import Sequence

class TestSequence(unittest.TestCase):

    def setUp(self):
        self.seq1 = Sequence("seq1", "GATTACA")

    def add_mock_gene(self, name="foo_gene"):
        mockgene = Mock()
        mockgene.identifier = name
        mockgene.indices = [2, 4]
        mockgene.death_flagged = False
        mockgene.to_mrna_fasta.return_value = "mockgene_to_mrna_fasta\n"
        mockgene.to_cds_fasta.return_value = "mockgene_to_cds_fasta\n"
        mockgene.to_protein_fasta.return_value = "mockgene_to_protein_fasta\n"
        mockgene.get_valid_mrnas = Mock(return_value=[])
        self.seq1.add_gene(mockgene)
        
    def add_mock_gene_with_1_mrna(self, name):
        mockgene = Mock()
        mockgene.indices = [1, 10]
        mockgene.identifier = name
        mockgene.death_flagged = False
        mockgene.mrnas = [Mock()]
        mockgene.mrnas[0].identifier = name+"-RA"
        mockgene.mrnas[0].cds = Mock()
        mockgene.mrnas[0].cds.identifier = [name+"-RA:CDS"]
        mockgene.mrnas[0].cds.length = Mock(return_value=5)
        mockgene.mrnas[0].exon = Mock()
        mockgene.mrnas[0].length = Mock(return_value=2)
        mockgene.get_valid_mrnas = Mock(return_value=mockgene.mrnas)
        mockgene.length = Mock(return_value=20)
        mockgene.get_partial_info.return_value = {"complete": 1, "start_no_stop": 0, "stop_no_start": 1, "no_stop_no_start": 1}
        mockgene.get_num_exons.return_value = 5
        mockgene.get_num_introns.return_value = 4
        mockgene.get_longest_exon.return_value = 20
        mockgene.get_longest_intron.return_value = 20
        mockgene.get_shortest_exon.return_value = 8
        mockgene.get_shortest_intron.return_value = 8
        mockgene.get_total_exon_length.return_value = 15
        mockgene.get_total_intron_length.return_value = 15
        self.seq1.add_gene(mockgene)
        
    def add_mock_gene_with_2_mrnas(self, name):
        mockgene = Mock()
        mockgene.indices = [20, 30]
        mockgene.identifier = name
        mockgene.death_flagged = False
        mockgene.mrnas = [Mock(), Mock()]
        mockgene.mrnas[0].identifier = name+"-RA"
        mockgene.mrnas[0].cds = None
        mockgene.mrnas[0].exon = None
        mockgene.mrnas[0].length = Mock(return_value=5)
        mockgene.mrnas[1].identifier = name+"-RB"
        mockgene.mrnas[1].cds = Mock()
        mockgene.mrnas[1].cds.identifier = [name+"-RB:CDS"]
        mockgene.mrnas[1].cds.length = Mock(return_value=3)
        mockgene.mrnas[1].exon = Mock()
        mockgene.mrnas[1].length = Mock(return_value=2)
        mockgene.get_valid_mrnas = Mock(return_value=mockgene.mrnas)
        mockgene.length = Mock(return_value=10)
        mockgene.get_partial_info.return_value = {"complete": 0, "start_no_stop": 1, "stop_no_start": 1, "no_stop_no_start": 1}
        mockgene.get_num_exons.return_value = 4
        mockgene.get_num_introns.return_value = 3
        mockgene.get_longest_exon.return_value = 10
        mockgene.get_longest_intron.return_value = 10
        mockgene.get_shortest_exon.return_value = 5
        mockgene.get_shortest_intron.return_value = 5
        mockgene.get_total_exon_length.return_value = 25
        mockgene.get_total_intron_length.return_value = 25
        self.seq1.add_gene(mockgene)

    def test_string(self):
        expected = "Sequence seq1 of length 7 containing 0 genes\n"
        self.assertEquals(expected, str(self.seq1))

    def test_how_many_Ns_forward(self):
        badseq = Sequence('seq1', 'NNnNNGATTACA')
        self.assertEqual(5, badseq.how_many_Ns_forward(1))

    def test_how_many_Ns_forward_returns_zero_if_no_Ns(self):
        badseq = Sequence('seq2', 'GATTACA')
        self.assertEqual(0, badseq.how_many_Ns_forward(3))

    def test_how_many_Ns_backward(self):
        badseq = Sequence('seq3', 'gattaNnN')
        self.assertEqual(3, badseq.how_many_Ns_backward(8))

    def test_how_many_Ns_backward_returns_zero_if_no_Ns(self):
        self.assertEqual(0, self.seq1.how_many_Ns_backward(3))

    def test_number_of_gagflags(self):
        gene1, gene2 = Mock(), Mock()
        gene1.number_of_gagflags.return_value = 2
        gene2.number_of_gagflags.return_value = 1
        self.seq1.genes = [gene1, gene2]
        self.assertEquals(3, self.seq1.number_of_gagflags())

    def test_remove_terminal_ns_beginning(self):
        badseq = Sequence('badseq', 'nNGATTACA')
        mockgene = Mock()
        mockgene.indices = [3, 6]
        badseq.genes = [mockgene]
        badseq.remove_terminal_ns()
        self.assertEquals("GATTACA", badseq.bases)

    def test_remove_terminal_ns_end(self):
        badseq = Sequence('badseq', 'GATTACAnNNn')
        mockgene = Mock()
        mockgene.indices = [2, 6]
        badseq.genes = [mockgene]
        badseq.remove_terminal_ns()
        self.assertEquals("GATTACA", badseq.bases)

    def test_remove_terminal_ns_beginning_and_end(self):
        badseq = Sequence('badseq', 'nnGATTACAnNNn')
        mockgene = Mock()
        mockgene.indices = [3, 8]
        badseq.genes = [mockgene]
        badseq.remove_terminal_ns()
        self.assertEquals("GATTACA", badseq.bases)

    def test_add_gene(self):
        self.add_mock_gene()
        self.assertEqual(1, len(self.seq1.genes))
    
    def test_remove_gene(self):
        self.add_mock_gene('foo_gene')
        self.assertEqual(1, len(self.seq1.genes))
        self.assertEqual(0, len(self.seq1.removed_genes))
        self.seq1.remove_gene('foo_gene')
        self.assertEqual(0, len(self.seq1.genes))
        self.assertEqual(1, len(self.seq1.removed_genes))

    def test_remove_genes_from_list(self):
        self.add_mock_gene('foo_gene')
        self.add_mock_gene('bar_gene')
        self.add_mock_gene('zub_gene')
        bad_genes = ["zub_gene", "foo_gene"]
        self.assertEquals(3, len(self.seq1.genes))
        removed_genes = self.seq1.remove_genes_from_list(bad_genes)
        self.assertEquals(2, len(removed_genes))
        self.assertEquals(1, len(self.seq1.genes))
        self.assertEquals(2, len(self.seq1.removed_genes))

    def test_remove_genes_from_list_bad_list(self):
        self.add_mock_gene('foo_gene')
        self.add_mock_gene('bar_gene')
        self.add_mock_gene('zub_gene')
        bad_genes = ["nice_gene", "bacon", 28]
        self.assertEquals(3, len(self.seq1.genes))
        self.seq1.remove_genes_from_list(bad_genes) # nothing should happen
        self.assertEquals(3, len(self.seq1.genes))

    def test_remove_mrnas_from_list(self):
        self.seq1.genes = [Mock()]
        self.seq1.genes[0].remove_mrnas_from_list.return_value = ["foo"]
        bad_mrnas = ["foo_mrna", "bar_mrna"]
        removed = self.seq1.remove_mrnas_from_list(bad_mrnas)
        self.assertEquals(["foo"], removed)
        self.seq1.genes[0].remove_mrnas_from_list.assert_called_with(bad_mrnas)
    
    def test_remove_empty_genes(self):
        self.add_mock_gene('foo_gene')
        self.add_mock_gene('bar_gene')
        self.add_mock_gene('zub_gene')
        self.seq1.genes[0].mrnas = [Mock()]
        self.seq1.genes[1].mrnas = []
        self.seq1.genes[2].mrnas = []
        self.assertEquals(3, len(self.seq1.genes))
        removed_genes = self.seq1.remove_empty_genes()
        self.assertEquals(2, len(removed_genes))
        self.assertEquals(1, len(self.seq1.genes))
        self.assertEquals(2, len(self.seq1.removed_genes))
    
    def test_remove_empty_mrnas(self):
        self.seq1.genes = [Mock(), Mock()]
        self.seq1.genes[0].remove_empty_mrnas.return_value = []
        self.seq1.genes[1].remove_empty_mrnas.return_value = []
        removed_mrnas = self.seq1.remove_empty_mrnas()
        self.seq1.genes[0].remove_empty_mrnas.assert_called_with()
        self.seq1.genes[1].remove_empty_mrnas.assert_called_with()

    def test_remove_empty_mrnas_returns_list(self):
        gene = Mock()
        gene.remove_empty_mrnas.return_value = [1, 2] #should be list of mRNAs but whatever
        self.seq1.genes = [gene]
        removed_mrnas = self.seq1.remove_empty_mrnas()
        self.assertEquals([1, 2], removed_mrnas)

    def test_remove_empty_mrnas_returns_list_multiple_genes(self):
        gene1 = Mock()
        gene1.remove_empty_mrnas.return_value = [1, 2]
        gene2 = Mock()
        gene2.remove_empty_mrnas.return_value = [3, 4]
        self.seq1.genes = [gene1, gene2]
        removed_mrnas = self.seq1.remove_empty_mrnas()
        self.assertEquals([1, 2, 3, 4], removed_mrnas)

    def test_get_gene_ids(self):
        self.seq1.genes = [Mock(), Mock()]
        self.seq1.genes[0].identifier = "foo gene"
        self.seq1.genes[1].identifier = "bar gene"
        expected = ["foo gene", "bar gene"]
        self.assertEquals(self.seq1.get_gene_ids(), expected)

    def test_get_mrna_ids(self):
        self.seq1.genes = [Mock(), Mock()]
        self.seq1.genes[0].get_mrna_ids.return_value = ["mrna1", "mrna2"]
        self.seq1.genes[1].get_mrna_ids.return_value = ["mrna3", "mrna4"]
        expected = ["mrna1", "mrna2", "mrna3", "mrna4"]
        self.assertEquals(self.seq1.get_mrna_ids(), expected)

    def test_trim_region(self):
        self.assertEquals("GATTACA", self.seq1.bases)
        self.seq1.trim_region(1, 4)
        self.assertEquals("ACA", self.seq1.bases)

    def test_trim_region_removes_gene_contained_in_trimmed_region(self):
        self.add_mock_gene()
        self.assertEquals(1, len(self.seq1.genes))
        self.seq1.trim_region(1, 3)
        self.assertEquals(0, len(self.seq1.genes))

    def test_add_annotations_from_list_adds_to_mrna(self):
        gene = Mock()
        mrna = Mock()
        self.seq1.genes = [gene]
        gene.mrnas = [mrna]
        gene.identifier = "foo_gene"
        gene.contains_mrna.return_value = True
        anno_list = [["foo_mrna", "Dbxref", "PFAM:0001"]]
        self.seq1.add_annotations_from_list(anno_list)
        gene.add_mrna_annotation.assert_called_with("foo_mrna", "Dbxref", "PFAM:0001")

    def test_add_annotations_from_list_adds_to_gene(self):
        gene = Mock()
        self.seq1.genes = [gene]
        gene.identifier = "foo_gene"
        anno_list = [["foo_gene", "name", "ABC123"], ["bar_gene", "name", "XYZ789"]]
        self.seq1.add_annotations_from_list(anno_list)
        self.assertEquals("ABC123", gene.name)

    def test_get_subseq(self):
        self.assertEquals("ATTA", self.seq1.get_subseq(2, 5))

    def test_get_cds_partial_info(self):
        self.add_mock_gene_with_1_mrna("foo_gene1")
        self.add_mock_gene_with_2_mrnas("foo_gene2")
        partial_info = self.seq1.get_cds_partial_info()
        self.assertEquals(1, partial_info["CDS: complete"])
    
    def test_get_contained_genes(self):
        fake_gene0 = Mock()
        fake_gene0.indices = [0, 10]
        self.seq1.add_gene(fake_gene0)
        fake_gene1 = Mock()
        fake_gene1.indices = [1, 9]
        self.seq1.add_gene(fake_gene1)
        fake_gene2 = Mock()
        fake_gene2.indices = [0, 10]
        self.seq1.add_gene(fake_gene2)
        fake_gene3 = Mock()
        fake_gene3.indices = [5, 15]
        self.seq1.add_gene(fake_gene3)
        fake_gene4 = Mock()
        fake_gene4.indices = [20, 30]
        self.seq1.add_gene(fake_gene4)
        contained = self.seq1.get_contained_genes()
        self.assertEqual(contained, [fake_gene1])
    
    def test_get_overlapping_genes(self):
        fake_gene0 = Mock()
        fake_gene0.indices = [0, 10]
        self.seq1.add_gene(fake_gene0)
        fake_gene1 = Mock()
        fake_gene1.indices = [1, 9]
        self.seq1.add_gene(fake_gene1)
        fake_gene2 = Mock()
        fake_gene2.indices = [0, 10]
        self.seq1.add_gene(fake_gene2)
        fake_gene3 = Mock()
        fake_gene3.indices = [5, 15]
        self.seq1.add_gene(fake_gene3)
        fake_gene4 = Mock()
        fake_gene4.indices = [20, 30]
        self.seq1.add_gene(fake_gene4)
        contained = self.seq1.get_overlapping_genes()
        self.assertTrue(fake_gene0 in contained)
        self.assertTrue(fake_gene1 in contained)
        self.assertTrue(fake_gene2 in contained)
        self.assertTrue(fake_gene3 in contained)
        

    def test_cds_to_gff(self):
        mockgene = Mock()
        mockgene.contains_mrna.return_value = True
        self.seq1.genes = [mockgene]
        foo = self.seq1.cds_to_gff("foo_mrna")
        mockgene.cds_to_gff.assert_called_with("seq1", "foo_mrna")

    def test_cds_to_tbl(self):
        mockgene = Mock()
        mockgene.contains_mrna.return_value = True
        self.seq1.genes = [mockgene]
        foo = self.seq1.cds_to_tbl("foo_mrna")
        mockgene.cds_to_tbl.assert_called_with("foo_mrna")

    def test_to_mrna_fasta(self):
        self.add_mock_gene()
        expected = "mockgene_to_mrna_fasta\n"
        self.assertEquals(expected, self.seq1.to_mrna_fasta())

    def test_to_cds_fasta(self):
        self.add_mock_gene()
        expected = "mockgene_to_cds_fasta\n"
        self.assertEquals(expected, self.seq1.to_cds_fasta())

    def test_to_protein_fasta(self):
        self.add_mock_gene()
        expected = "mockgene_to_protein_fasta\n"
        self.assertEquals(expected, self.seq1.to_protein_fasta())

    def test_to_tbl(self):
        self.add_mock_gene()
        self.seq1.genes[0].to_tbl.return_value = "mockgene to tbl"
        tbl = self.seq1.to_tbl()
        expected = ">Feature seq1\n"
        expected += "1\t7\tREFERENCE\n"
        expected += "\t\t\tPBARC\t12345\n"
        expected += "mockgene to tbl"
        self.assertEquals(tbl, expected)

    def test_stats(self):
        self.add_mock_gene_with_1_mrna("foo_gene1")
        self.add_mock_gene_with_2_mrnas("foo_gene2")
        stats = self.seq1.stats()
        self.assertEquals(stats["Total sequence length"], 7)
        self.assertEquals(stats["Number of genes"], 2)
        self.assertEquals(stats["Number of mRNAs"], 3)
        self.assertEquals(stats["Number of exons"], 9)
        self.assertEquals(stats["Number of introns"], 7)
        self.assertEquals(stats["Number of CDS"], 2)
        self.assertEquals(stats["Overlapping genes"], 0)
        self.assertEquals(stats["Contained genes"], 0)
        self.assertEquals(stats["CDS: complete"], 1)
        self.assertEquals(stats["CDS: start, no stop"], 1)
        self.assertEquals(stats["CDS: stop, no start"], 2)
        self.assertEquals(stats["CDS: no stop, no start"], 2)
        self.assertEquals(stats["Longest gene"], 20)
        self.assertEquals(stats["Longest mRNA"], 5)
        self.assertEquals(stats["Longest exon"], 20)
        self.assertEquals(stats["Longest intron"], 20)
        self.assertEquals(stats["Longest CDS"], 5)
        self.assertEquals(stats["Shortest gene"], 10)
        self.assertEquals(stats["Shortest mRNA"], 2)
        self.assertEquals(stats["Shortest exon"], 5)
        self.assertEquals(stats["Shortest intron"], 5)
        self.assertEquals(stats["Shortest CDS"], 3)
        self.assertEquals(stats["Total gene length"], 30)
        self.assertEquals(stats["Total mRNA length"], 9)
        self.assertEquals(stats["Total exon length"], 40)
        self.assertEquals(stats["Total intron length"], 40)
        self.assertEquals(stats["Total CDS length"], 8)



##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestSequence))
    return suite

if __name__ == '__main__':
    unittest.main()
