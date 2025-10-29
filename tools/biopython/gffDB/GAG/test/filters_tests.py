#!/usr/bin/env python

import unittest
from mock import Mock
from src.gene import Gene
from src.sequence import Sequence
from src.filters import *

class TestFilters(unittest.TestCase):
        
    def test_min_cds_length_filter(self):
        cds_length = MinCDSLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        
        # Give the sequence some genes
        seq.genes = [Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo1'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo2'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo3')]
        
        # Give the mock mrnas some cds's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.cds = Mock()
        test_mrna0.cds.length = Mock(return_value=20)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.cds = None
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.cds = Mock()
        test_mrna2.cds.length = Mock(return_value=30)
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.cds = Mock()
        test_mrna3.cds.length = Mock(return_value=40)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[0].death_flagged = False
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[1].death_flagged = False
        seq.genes[2].mrnas = [test_mrna3]
        seq.genes[2].death_flagged = False
        
        # Apply the filter
        cds_length.apply(seq)

        self.assertEqual(len(seq.genes), 3)
        self.assertEqual(seq.genes[1].mrnas, [test_mrna1, test_mrna2])
        self.assertEqual(seq.genes[2].mrnas, [test_mrna3])

###################################################################################################
        
    def test_max_cds_length_filter(self):
        cds_length = MaxCDSLengthFilter(100)
    
        # Create a mock sequence
        seq = Sequence()
        
        # Give the sequence some genes
        seq.genes = [Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo1'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo2'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo3')]
        
        # Give the mock mrnas some cds's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.cds = Mock()
        test_mrna0.cds.length = Mock(return_value=90)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.cds = None
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.cds = Mock()
        test_mrna2.cds.length = Mock(return_value=100)
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.cds = Mock()
        test_mrna3.cds.length = Mock(return_value=110)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[0].death_flagged = False
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[1].death_flagged = False
        seq.genes[2].mrnas = [test_mrna3]
        seq.genes[2].death_flagged = False
        
        # Apply the filter
        cds_length.apply(seq)

        self.assertEqual(len(seq.genes), 3)
        self.assertEqual(seq.genes[0].mrnas, [test_mrna0])
        self.assertEqual(seq.genes[1].mrnas, [test_mrna1, test_mrna2])

###################################################################################################
        
    def test_min_exon_length_filter(self):
        exon_length = MinExonLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        
        # Give the sequence some genes
        seq.genes = [Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo1'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo2'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo3')]
        
        # Give the mock mrnas some exon's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.exon = Mock()
        test_mrna0.get_shortest_exon = Mock(return_value=20)
        test_mrna0.get_longest_exon = Mock(return_value=20)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.exon = Mock()
        test_mrna1.get_shortest_exon = Mock(return_value=30)
        test_mrna1.get_longest_exon = Mock(return_value=30)
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.exon = None
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.exon = Mock()
        test_mrna3.get_shortest_exon = Mock(return_value=40)
        test_mrna3.get_longest_exon = Mock(return_value=40)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[0].death_flagged = False
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[1].death_flagged = False
        seq.genes[2].mrnas = [test_mrna3]
        seq.genes[2].death_flagged = False
        
        # Apply the filter
        exon_length.apply(seq)
        
        self.assertEqual(len(seq.genes), 3)
        self.assertEqual(seq.genes[1].mrnas, [test_mrna1, test_mrna2])
        self.assertEqual(seq.genes[2].mrnas, [test_mrna3])

###################################################################################################

    def test_max_exon_length_filter(self):
        exon_length = MaxExonLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        
        # Give the sequence some genes
        seq.genes = [Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo1'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo2'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo3')]
        
        # Give the mock mrnas some exon's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.exon = Mock()
        test_mrna0.get_shortest_exon = Mock(return_value=20)
        test_mrna0.get_longest_exon = Mock(return_value=20)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.exon = Mock()
        test_mrna1.get_shortest_exon = Mock(return_value=30)
        test_mrna1.get_longest_exon = Mock(return_value=30)
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.exon = None
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.exon = Mock()
        test_mrna3.get_shortest_exon = Mock(return_value=40)
        test_mrna3.get_longest_exon = Mock(return_value=40)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[0].death_flagged = False
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[1].death_flagged = False
        seq.genes[2].mrnas = [test_mrna3]
        seq.genes[2].death_flagged = False
        
        # Apply the filter
        exon_length.apply(seq)
        
        self.assertEqual(len(seq.genes), 3)
        self.assertEqual(seq.genes[0].mrnas, [test_mrna0])
        self.assertEqual(seq.genes[1].mrnas, [test_mrna1, test_mrna2])

###################################################################################################

    def test_min_intron_length_filter(self):
        intron_length = MinIntronLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        
        # Give the sequence some genes
        seq.genes = [Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo1'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo2'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo3')]
        
        # Give the mock mrnas some exon's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.exon = Mock()
        test_mrna0.get_shortest_intron = Mock(return_value=20)
        test_mrna0.get_longest_intron = Mock(return_value=20)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.exon = Mock()
        test_mrna1.get_shortest_intron = Mock(return_value=30)
        test_mrna1.get_longest_intron = Mock(return_value=30)
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.exon = None
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.exon = Mock()
        test_mrna3.get_shortest_intron = Mock(return_value=40)
        test_mrna3.get_longest_intron = Mock(return_value=40)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[0].death_flagged = False
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[1].death_flagged = False
        seq.genes[2].mrnas = [test_mrna3]
        seq.genes[2].death_flagged = False
        
        # Apply the filter
        intron_length.apply(seq)
        
        self.assertEqual(len(seq.genes), 3)
        self.assertEqual(seq.genes[1].mrnas, [test_mrna1, test_mrna2])
        self.assertEqual(seq.genes[2].mrnas, [test_mrna3])

    def test_min_intron_length_filter_doesnt_remove_single_exon_mrnas(self):
        intron_length = MinIntronLengthFilter(30)
        seq = Sequence()
        gene, mrna, exon = Mock(), Mock(), Mock()
        gene.death_flagged = False
        seq.genes = [gene]
        gene.mrnas = [mrna]
        mrna.identifier = "badmrna"
        mrna.get_shortest_intron.return_value = 0
        mrna.death_flagged = False
        mrna.exon = exon
        self.assertEquals(1, len(seq.genes[0].mrnas))
        intron_length.apply(seq)
        # Filter shouldn't remove mRNA based on a shortest_intron = 0 return value
        self.assertEquals(1, len(seq.genes[0].mrnas))
        assert not gene.remove_mrna.called

###################################################################################################

    def test_max_intron_length_filter(self):
        intron_length = MaxIntronLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        
        # Give the sequence some genes
        seq.genes = [Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo1'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo2'), Gene('foo_seq', 'geib_labs', [1, 2], '+', 'foo3')]
        
        # Give the mock mrnas some exon's
        test_mrna0 = Mock()
        test_mrna0.identifier = 'foo1-RA'
        test_mrna0.death_flagged = False
        test_mrna0.exon = Mock()
        test_mrna0.get_shortest_intron = Mock(return_value=20)
        test_mrna0.get_longest_intron = Mock(return_value=20)
        
        test_mrna1 = Mock()
        test_mrna1.identifier = 'foo2-RA'
        test_mrna1.death_flagged = False
        test_mrna1.exon = Mock()
        test_mrna1.get_shortest_intron = Mock(return_value=30)
        test_mrna1.get_longest_intron = Mock(return_value=30)
        
        test_mrna2 = Mock()
        test_mrna2.identifier = 'foo2-RB'
        test_mrna2.death_flagged = False
        test_mrna2.exon = None
        
        test_mrna3 = Mock()
        test_mrna3.identifier = 'foo3-RA'
        test_mrna3.death_flagged = False
        test_mrna3.exon = Mock()
        test_mrna3.get_shortest_intron = Mock(return_value=40)
        test_mrna3.get_longest_intron = Mock(return_value=40)
        
        # Give the mock genes some mrnas
        seq.genes[0].mrnas = [test_mrna0]
        seq.genes[0].death_flagged = False
        seq.genes[1].mrnas = [test_mrna1, test_mrna2]
        seq.genes[1].death_flagged = False
        seq.genes[2].mrnas = [test_mrna3]
        seq.genes[2].death_flagged = False
        
        # Apply the filter
        intron_length.apply(seq)
        
        self.assertEqual(len(seq.genes), 3)
        self.assertEqual(seq.genes[0].mrnas, [test_mrna0])
        self.assertEqual(seq.genes[1].mrnas, [test_mrna1, test_mrna2])

###################################################################################################

    def test_min_gene_length_filter(self):
        gene_length_range = MinGeneLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        test_gene0 = Mock(identifier = 'foo1')
        test_gene1 = Mock(identifier = 'foo2')
        test_gene2 = Mock(identifier = 'foo3')
        
        test_gene0.death_flagged = False
        test_gene1.death_flagged = False
        test_gene2.death_flagged = False
        
        test_gene0.length = Mock(return_value=20)
        test_gene1.length = Mock(return_value=30)
        test_gene2.length = Mock(return_value=40)
        
        # Give the sequence some genes
        seq.genes = [test_gene0, test_gene1, test_gene2]
        
        # Apply the filter
        gene_length_range.apply(seq)
        
        self.assertEqual(seq.genes, [test_gene1, test_gene2])

###################################################################################################

    def test_max_gene_length_filter(self):
        gene_length_range = MaxGeneLengthFilter(30)
    
        # Create a mock sequence
        seq = Sequence()
        
        test_gene0 = Mock(identifier = 'foo1')
        test_gene1 = Mock(identifier = 'foo2')
        test_gene2 = Mock(identifier = 'foo3')
        
        test_gene0.death_flagged = False
        test_gene1.death_flagged = False
        test_gene2.death_flagged = False
        
        test_gene0.length = Mock(return_value=20)
        test_gene1.length = Mock(return_value=30)
        test_gene2.length = Mock(return_value=40)
        
        # Give the sequence some genes
        seq.genes = [test_gene0, test_gene1, test_gene2]
        
        # Apply the filter
        gene_length_range.apply(seq)
        
        self.assertEqual(seq.genes, [test_gene0, test_gene1])




##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFilters))
    return suite

if __name__ == '__main__':
    unittest.main()
