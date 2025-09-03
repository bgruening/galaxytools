#!/usr/bin/env python

import copy
from src.gene_part import *

class Exon(GenePart):

    def __init__(self, **kwargs):
        kwargs['feature_type'] = 'exon'
        GenePart.__init__(self, **kwargs)

    def to_tbl(self, has_start, has_stop, feature_type):
        """Returns a string representing the .tbl-formatted entry for this exon."""
        indices = copy.deepcopy(self.indices)
        return write_tbl_entry(indices, self.strand, has_start, has_stop, feature_type, self.annotations)

