#!/usr/bin/env python

import ast

# All the filters
from src.filters import *

class FilterManager:

    def __init__(self):
        # Build filters
        self.filters = dict()
        self.filters['cds_shorter_than'] = MinCDSLengthFilter()
        self.filters['cds_longer_than'] = MaxCDSLengthFilter()
        self.filters['exon_shorter_than'] = MinExonLengthFilter()
        self.filters['exon_longer_than'] = MaxExonLengthFilter()
        self.filters['intron_shorter_than'] = MinIntronLengthFilter()
        self.filters['intron_longer_than'] = MaxIntronLengthFilter()
        self.filters['gene_shorter_than'] = MinGeneLengthFilter()
        self.filters['gene_longer_than'] = MaxGeneLengthFilter()
    
    def apply_filter(self, filter_name, val, filter_mode, seq):
        val = ast.literal_eval(val)
        self.filters[filter_name].arg = val
        self.filters[filter_name].filter_mode = filter_mode 
        self.filters[filter_name].apply(seq)
    
    def get_filter_arg(self, filter_name):
        return self.filters[filter_name].arg
   
