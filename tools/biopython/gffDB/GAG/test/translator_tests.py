#!/usr/bin/env python

import unittest
from src.translator import *

class TestTranslate(unittest.TestCase):

    def test_translate_one_codon(self):
        self.assertEquals('T', translate('act', '+'))
        self.assertEquals('S', translate('act', '-'))

    def test_translate_longer_sequence(self):
        test_seq = 'CATGACAGAAGATATTTC'
        self.assertEquals('HDRRYF', translate(test_seq, '+',))
        self.assertEquals('EISSVM', translate(test_seq, '-'))
        
    def test_valid_seq(self):
        self.assertTrue(valid_seq('actg'))
        self.assertFalse(valid_seq('acth'))
        self.assertFalse(valid_seq('ac'))

    def test_has_start_codon(self):
        self.assertTrue(has_start_codon('auggattaca'))
        self.assertFalse(has_start_codon('guggattaca')) # currently no support for alternate start codons

    def test_has_stop_codon(self):
        self.assertTrue(has_stop_codon('gattacatag'))
        self.assertTrue(has_stop_codon('gattacataa'))
        self.assertTrue(has_stop_codon('gattacatga'))
        self.assertFalse(has_stop_codon('gattacaact'))

    def test_reverse_complement(self):
        self.assertEquals('C', reverse_complement('G'))
        self.assertEquals('CAT', reverse_complement('ATG'))

    def test_reverse_complement_with_bogus_base(self):
        self.assertEquals('CATN', reverse_complement('MATG'))

    def test_reverse_complement_longer_seq(self):
        self.assertEquals('TGTAATCTGTAATCTGTAATCTGTAATCTGTAATC', reverse_complement('GATTACAGATTACAGATTACAGATTACAGATTACA'))

    def test_translate_with_n_in_seq(self):
        test_seq = 'CATGACAGAAGATNTTTC'
        self.assertEquals('HDRRXF', translate(test_seq, '+'))

    def test_contains_internal_stop(self):
        test_seq = 'gattaggat' # translates to 'D*D'
        self.assertTrue(contains_internal_stop(test_seq, '+'))

    def test_contains_internal_stop_false(self):
        test_seq = 'GATTACTAG' # stop, but not internal
        self.assertFalse(contains_internal_stop(test_seq, '+'))

        
        
##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestTranslate))
    return suite

if __name__ == '__main__':
    unittest.main()
