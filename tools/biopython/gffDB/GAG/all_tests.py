#!/usr/bin/env python

# import all the lovely files
import unittest
import test.fasta_reader_tests
import test.gene_part_tests
import test.xrna_tests
import test.gene_tests
import test.translator_tests
import test.gff_reader_tests
import test.sequence_tests
import test.filter_manager_tests
import test.filters_tests
import test.stats_manager_tests
import test.seq_helper_tests
import test.cds_tests
import test.exon_tests

# get suites from test modules
suite1 = test.fasta_reader_tests.suite()
suite2 = test.gene_part_tests.suite()
suite3 = test.xrna_tests.suite()
suite4 = test.gene_tests.suite()
suite9 = test.translator_tests.suite()
suite10 = test.gff_reader_tests.suite()
suite11 = test.sequence_tests.suite()
suite12 = test.filter_manager_tests.suite()
suite13 = test.filters_tests.suite()
suite14 = test.stats_manager_tests.suite()
suite15 = test.seq_helper_tests.suite()
suite16 = test.cds_tests.suite()
suite17 = test.exon_tests.suite()

# collect suites in a TestSuite object
suite = unittest.TestSuite()
suite.addTest(suite1)
suite.addTest(suite2)
suite.addTest(suite3)
suite.addTest(suite4)
suite.addTest(suite9)
suite.addTest(suite10)
suite.addTest(suite11)
suite.addTest(suite12)
suite.addTest(suite13)
suite.addTest(suite14)
suite.addTest(suite15)
suite.addTest(suite16)
suite.addTest(suite17)

# run suite
unittest.TextTestRunner(verbosity=2).run(suite)
