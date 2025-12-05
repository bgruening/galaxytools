#!/usr/bin/env python

import unittest
import io
import os
from mock import Mock, patch, PropertyMock
from src.gff_reader import *

class TestGFFReader(unittest.TestCase):

    def setUp(self):
        self.reader = GFFReader()

    def tearDown(self):
        # Remove extra files created by GFFReader
        try:
            os.remove("genome.comments.gff")
            os.remove("genome.invalid.gff")
            os.remove("genome.ignored.gff")
        except OSError:
            pass

    def test_validate_line_not_enough_fields(self):
        badline = "scaffold00080\tmaker\tgene\t106151\t109853\t+\t.\tID=BDOR_007864\n"
        self.assertFalse(self.reader.validate_line(badline))

    def test_validate_line_no_id(self):
        badline = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tName=BDOR_007864\n"
        self.assertFalse(self.reader.validate_line(badline))

    def test_validate_line_indices_out_of_order(self):
        badline = "scaffold00080\tmaker\tgene\t109853\t106151\t.\t+\t.\tID=BDOR_007864;Name=BDOR_007864\n"
        self.assertFalse(self.reader.validate_line(badline))

    def test_validate_line(self):
        goodline = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n"
        self.assertTrue(self.reader.validate_line(goodline))

    def test_line_type_gene(self):
        line = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n".split('\t')
        self.assertEqual('gene', self.reader.line_type(line))

    def test_line_type_mrna(self):
        line = "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n".split('\t')
        self.assertEqual('mRNA', self.reader.line_type(line))

    def test_line_type_exon(self):
        line = "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('exon', self.reader.line_type(line))

    def test_line_type_cds(self):
        line = "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('CDS', self.reader.line_type(line))

    def test_line_type_start_codon(self):
        line = "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('start_codon', self.reader.line_type(line))

    def test_line_type_stop_codon(self):
        line = "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n".split('\t')
        self.assertEqual('stop_codon', self.reader.line_type(line))

    def test_parse_attributes(self):
        attr = "ID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n"
        parsed = self.reader.parse_attributes(attr)
        self.assertEqual('BDOR_007864-RA:stop2', parsed['identifier'])
        self.assertEqual('BDOR_007864-RA', parsed['parent_id'])

    def test_parse_attributes_with_name(self):
        attr = "ID=BDOR_007864-RA:stop2;Name=BDOR_007864-RA;Parent=BDOR_007864-RA\n"
        parsed = self.reader.parse_attributes(attr)
        self.assertEqual('BDOR_007864-RA:stop2', parsed['identifier'])
        self.assertEqual('BDOR_007864-RA', parsed['parent_id'])
        self.assertEqual('BDOR_007864-RA', parsed['name'])

    def test_extract_cds_args(self):
        line = "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n".split('\t')
        args = self.reader.extract_cds_args(line)
        expected = {'indices': [106151, 106451], 'strand': '+', 'phase': 0, 'identifier': 'BDOR_007864-RA:cds:0', 'parent_id': 'BDOR_007864-RA'}
        self.assertEqual(expected, args)

    def test_extract_exon_args(self):
        line = "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n".split('\t')
        expected = {'indices': [106151, 106451], 'score': 0.9, 'strand': '+', 'identifier': 'BDOR_007864-RA:exon:0', 'parent_id': 'BDOR_007864-RA'}
        args = self.reader.extract_exon_args(line)
        self.assertEqual(expected, args)

    def test_extract_mrna_args(self):
        line = "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n".split('\t')
        expected = {'indices': [106151, 109853], 'identifier': 'BDOR_007864-RA', 'strand': '+', 'parent_id': 'BDOR_007864',
                    'seq_name': "scaffold00080", 'source': "maker"}
        args = self.reader.extract_mrna_args(line)
        self.assertEqual(expected, args)

    def test_extract_gene_args(self):
        line = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n".split('\t')
        expected = {'seq_name': 'scaffold00080', 'source': 'maker', 'indices': [106151, 109853],\
                    'strand': '+', 'identifier': 'BDOR_007864'}
        args = self.reader.extract_gene_args(line)
        self.assertEqual(expected, args)

    def test_extract_other_feature_args(self):
        line = "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n".split('\t')
        expected = {'feature_type': 'start_codon', 'indices': [106151, 106153],\
                    'identifier': 'BDOR_007864-RA:start1', 'parent_id': 'BDOR_007864-RA'}
        args = self.reader.extract_other_feature_args(line)
        self.assertEqual(expected, args)

    def test_update_cds(self):
        current_cds = Mock()
        line = "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2\tID=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RA\n".split('\t')
        self.reader.update_cds(line, current_cds)
        current_cds.add_indices.assert_called_with([106509, 106749])
        current_cds.add_phase.assert_called_with(2)
        current_cds.add_identifier.assert_called_with('BDOR_007864-RA:cds:1')

    def test_update_exon(self):
        current_exon = Mock()
        line = "scaffold00080\tmaker\texon\t106509\t106749\t8.34\t+\t2\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RA\n".split('\t')
        self.reader.update_exon(line, current_exon)
        current_exon.add_indices.assert_called_with([106509, 106749])
        current_exon.add_identifier.assert_called_with('BDOR_007864-RA:exon:1')

    def test_read_file(self):
        text = self.get_sample_text()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertEquals(2, len(genes))
        self.assertEquals('BDOR_007864-RA', genes[0].mrnas[0].identifier)
        self.assertEquals([179489, 179691], genes[1].mrnas[0].cds.indices[2])

    def get_sample_text(self):
        sample_text = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\texon\t106509\t106749\t0.9\t+\t.\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2\tID=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tgene\t145206\t183302\t.\t+\t.\tID=BDOR_007866\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t145206\t183302\t.\t+\t.\tID=BDOR_007866-RB;Parent=BDOR_007866\n"
        sample_text += "scaffold00080\tmaker\texon\t145206\t145282\t0.065333\t+\t.\tID=BDOR_007866-RB:exon:5;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\texon\t145607\t145865\t47.919\t+\t.\tID=BDOR_007866-RB:exon:6;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\texon\t145928\t146176\t67.378\t+\t.\tID=BDOR_007866-RB:exon:7;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t145206\t145282\t.\t+\t.\tID=BDOR_007866-RB:UTR1;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t145607\t145865\t.\t+\t.\tID=BDOR_007866-RB:UTR2;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t145928\t146176\t.\t+\t.\tID=BDOR_007866-RB:UTR3;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tfive_prime_UTR\t154498\t154575\t.\t+\t.\tID=BDOR_007866-RB:UTR4;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t154576\t154620\t.\t+\t0\tID=BDOR_007866-RB:cds:5;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t179210\t179419\t.\t+\t0\tID=BDOR_007866-RB:cds:6;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t179489\t179691\t.\t+\t0\tID=BDOR_007866-RB:cds:7;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tthree_prime_UTR\t183025\t183302\t.\t+\t.\tID=BDOR_007866-RB:UTR5;Parent=BDOR_007866-RB\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t183022\t183024\t.\t+\t.\tID=BDOR_007866-RB:stop2;Parent=BDOR_007866-RB\n"
        return sample_text

    def get_out_of_order_text(self):
        sample_text = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RB;Parent=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\texon\t106509\t106749\t0.9\t+\t.\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2\tID=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\texon\t106509\t106749\t0.9\t+\t.\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2\tID=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RB\n"
        return sample_text
        
    def get_out_of_order_text_with_missing_parent(self):
        sample_text = "scaffold00080\tmaker\tgene\t106151\t109853\t.\t+\t.\tID=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\tmRNA\t106151\t109853\t.\t+\t.\tID=BDOR_007864-RA;Parent=BDOR_007864\n"
        sample_text += "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\texon\t106509\t106749\t0.9\t+\t.\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2\tID=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RA\n"
        sample_text += "scaffold00080\tmaker\texon\t106151\t106451\t0.9\t+\t.\tID=BDOR_007864-RA:exon:0;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\texon\t106509\t106749\t0.9\t+\t.\tID=BDOR_007864-RA:exon:1;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106151\t106451\t.\t+\t0\tID=BDOR_007864-RA:cds:0;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tCDS\t106509\t106749\t.\t+\t2\tID=BDOR_007864-RA:cds:1;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tstart_codon\t106151\t106153\t.\t+\t.\tID=BDOR_007864-RA:start1;Parent=BDOR_007864-RB\n"
        sample_text += "scaffold00080\tmaker\tstop_codon\t109851\t109853\t.\t+\t.\tID=BDOR_007864-RA:stop2;Parent=BDOR_007864-RB\n"
        return sample_text
        
    def test_read_file_out_of_order(self):
        text = self.get_out_of_order_text()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertEqual(1, len(genes))
        self.assertEqual('BDOR_007864-RA', genes[0].mrnas[0].identifier)
        self.assertEqual(2, len(genes[0].mrnas))
        self.assertEqual(2, len(genes[0].mrnas[0].exon.indices))
        self.assertEqual(2, len(genes[0].mrnas[1].exon.indices))

    def test_read_file_doesnt_loop_infinitely_when_feature_with_no_parent_mrna(self):
        text = self.get_out_of_order_text_with_missing_parent()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertEqual(1, len(genes))
        
    def get_annotated_gff(self):
        result = "Scaffold1\tI5K\tgene\t133721\t162851\t.\t-\t.\tID=AGLA000002;Name=AglaTmpM000002;\n"
        result += "Scaffold1\tI5K\tmRNA\t133721\t162851\t.\t-\t.\tID=AGLA000002-RA;Name=AglaTmpM000002-RA;Parent=AGLA000002;Dbxref=PRINTS:PR00075;\n"
        result += "Scaffold1\tI5K\texon\t133721\t135519\t.\t-\t.\tID=AGLA000002-RA-EXON01;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t140163\t140635\t.\t-\t.\tID=AGLA000002-RA-EXON02;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t147266\t147396\t.\t-\t.\tID=AGLA000002-RA-EXON03;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t152757\t152979\t.\t-\t.\tID=AGLA000002-RA-EXON04;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t162720\t162762\t.\t-\t.\tID=AGLA000002-RA-EXON05;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t162825\t162851\t.\t-\t.\tID=AGLA000002-RA-EXON06;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t140426\t140635\t.\t-\t0\tID=AGLA000002-RA-CDS01;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t147266\t147396\t.\t-\t2\tID=AGLA000002-RA-CDS02;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t152757\t152976\t.\t-\t0\tID=AGLA000002-RA-CDS03;Parent=AGLA000002-RA;\n"
        return result

    def get_annotated_gff_multi_dbxref(self):
        result = "Scaffold1\tI5K\tgene\t133721\t162851\t.\t-\t.\tID=AGLA000002;Name=AglaTmpM000002;\n"
        result += "Scaffold1\tI5K\tmRNA\t133721\t162851\t.\t-\t.\tID=AGLA000002-RA;Name=AglaTmpM000002-RA;Parent=AGLA000002;Dbxref=PRINTS:PR00075,PFAM:foo;\n"
        result += "Scaffold1\tI5K\texon\t133721\t135519\t.\t-\t.\tID=AGLA000002-RA-EXON01;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t140163\t140635\t.\t-\t.\tID=AGLA000002-RA-EXON02;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t147266\t147396\t.\t-\t.\tID=AGLA000002-RA-EXON03;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t152757\t152979\t.\t-\t.\tID=AGLA000002-RA-EXON04;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t162720\t162762\t.\t-\t.\tID=AGLA000002-RA-EXON05;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t162825\t162851\t.\t-\t.\tID=AGLA000002-RA-EXON06;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t140426\t140635\t.\t-\t0\tID=AGLA000002-RA-CDS01;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t147266\t147396\t.\t-\t2\tID=AGLA000002-RA-CDS02;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t152757\t152976\t.\t-\t0\tID=AGLA000002-RA-CDS03;Parent=AGLA000002-RA;\n"
        return result

    def get_annotated_gff_multi_dbxref_repeated_anno(self):
        result = "Scaffold1\tI5K\tgene\t133721\t162851\t.\t-\t.\tID=AGLA000002;Name=AglaTmpM000002;\n"
        result += "Scaffold1\tI5K\tmRNA\t133721\t162851\t.\t-\t.\tID=AGLA000002-RA;Name=AglaTmpM000002-RA;Parent=AGLA000002;Dbxref=PRINTS:PR00075;Dbxref=PFAM:foo;\n"
        result += "Scaffold1\tI5K\texon\t133721\t135519\t.\t-\t.\tID=AGLA000002-RA-EXON01;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t140163\t140635\t.\t-\t.\tID=AGLA000002-RA-EXON02;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t147266\t147396\t.\t-\t.\tID=AGLA000002-RA-EXON03;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t152757\t152979\t.\t-\t.\tID=AGLA000002-RA-EXON04;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t162720\t162762\t.\t-\t.\tID=AGLA000002-RA-EXON05;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\texon\t162825\t162851\t.\t-\t.\tID=AGLA000002-RA-EXON06;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t140426\t140635\t.\t-\t0\tID=AGLA000002-RA-CDS01;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t147266\t147396\t.\t-\t2\tID=AGLA000002-RA-CDS02;Parent=AGLA000002-RA;\n"
        result += "Scaffold1\tI5K\tCDS\t152757\t152976\t.\t-\t0\tID=AGLA000002-RA-CDS03;Parent=AGLA000002-RA;\n"
        return result

    def test_read_file_annotated(self):
        text = self.get_annotated_gff()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertEquals(1, len(genes))
        self.assertEquals({"Dbxref": ["PRINTS:PR00075"]}, genes[0].mrnas[0].annotations)

    def test_read_file_annotated_multi_dbxref(self):
        text = self.get_annotated_gff_multi_dbxref()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertEquals(1, len(genes))
        self.assertEquals({"Dbxref": ["PRINTS:PR00075", "PFAM:foo"]}, genes[0].mrnas[0].annotations)

    def test_read_file_annotated_multi_dbxref_repeated_anno(self):
        text = self.get_annotated_gff_multi_dbxref_repeated_anno()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertEquals(1, len(genes))
        self.assertEquals({"Dbxref": ["PRINTS:PR00075", "PFAM:foo"]}, genes[0].mrnas[0].annotations)

    def test_CDS_knows_its_strand(self):
        text = self.get_annotated_gff()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertTrue(genes[0].mrnas[0].cds)
        self.assertEquals('-', genes[0].mrnas[0].cds.strand)

    def test_exon_knows_its_strand(self):
        text = self.get_annotated_gff()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertTrue(genes[0].mrnas[0].exon)
        self.assertEquals('-', genes[0].mrnas[0].exon.strand)

    def test_mrna_knows_its_strand(self):
        text = self.get_annotated_gff()
        inbuff = io.BytesIO(text)
        genes, comments, invalids, ignored = self.reader.read_file(inbuff)
        self.assertTrue(genes[0].mrnas[0])
        self.assertEquals('-', genes[0].mrnas[0].strand)


##########################
def suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(TestGFFReader))
    return suite

if __name__ == '__main__':
    unittest.main()
