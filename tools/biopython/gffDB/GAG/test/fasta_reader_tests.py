#!/usr/bin/env python

import unittest
import io
from src.fasta_reader import FastaReader
from mock import Mock

class TestFastaReader(unittest.TestCase):

    def setUp(self):
        self.reader = FastaReader()

    def test_read(self):
        no_line_breaks = io.BytesIO('>seq_1\nGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA\n' +
                                    '>seq_2\nNNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN')
        line_breaks = io.BytesIO('>seq_1\nGATTACAGATTACAGATTACAGATTACA\nGATTACAGATTACAGATTACAGATTACA\n' +
                                 '>seq_2\nNNNNNNNNGATTACAGATTACAGATTAC\nANNNNNNNNNNN')

        self.reader.read(no_line_breaks)
        self.assertEquals(2, len(self.reader.seqs))
        self.assertEquals('seq_1', self.reader.seqs[0].header)
        self.assertEquals('GATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACAGATTACA', self.reader.seqs[0].bases)
        self.assertEquals('seq_2', self.reader.seqs[1].header)
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', self.reader.seqs[1].bases)
        self.reader.read(line_breaks)
        self.assertEquals(4, len(self.reader.seqs))
        self.assertEquals('NNNNNNNNGATTACAGATTACAGATTACANNNNNNNNNNN', self.reader.seqs[3].bases)

        


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestFastaReader))
    return suite

if __name__ == '__main__':
    unittest.main()
