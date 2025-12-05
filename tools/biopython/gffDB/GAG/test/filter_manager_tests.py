#!/usr/bin/env python

import unittest
from mock import Mock
from src.filter_manager import FilterManager

class TestFilterManager(unittest.TestCase):

    def setUp(self):
        self.filter_mgr = FilterManager()
        
    def test_modify_filter_arg(self):
        fake_seq = Mock()
        self.filter_mgr.filters['cds_shorter_than'] = Mock()
        self.filter_mgr.filters['cds_shorter_than'].arg = 0
        self.filter_mgr.filters['cds_shorter_than'].remove = True
        self.filter_mgr.apply_filter('cds_shorter_than', '30', False, fake_seq)
        self.assertEqual(self.filter_mgr.get_filter_arg('cds_shorter_than'), 30)
        self.filter_mgr.filters['cds_shorter_than'].apply.assert_called_with(fake_seq)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFilterManager))
    return suite

if __name__ == '__main__':
    unittest.main()
