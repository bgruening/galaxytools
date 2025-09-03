#!/usr/bin/env python

import copy
from src.gene_part import *

class CDS(GenePart):

    def __init__(self, identifier=None, indices=None, \
                 score=None, phase=None, strand=None, parent_id=None):
        GenePart.__init__(self, feature_type='CDS', identifier=identifier, \
                indices=indices, score=score, strand=strand, parent_id=parent_id)
        self.phase = []
        if phase is not None:
            self.phase.append(phase)

    def get_phase(self, i):
        """Returns phase for given segment of CDS."""
        if self.phase and len(self.phase) > i:
            return self.phase[i]
        else:
            return "."

    def add_phase(self, ph):
        """Appends phase to CDS"""
        self.phase.append(ph)

    def get_start_indices(self, strand):
        """Returns coordinates of first and third base of CDS."""
        if strand == '+':
            first_index = self.indices[0][0]
            return [first_index, first_index+2]
        elif strand == '-':
            first_index = self.indices[0][1]
            return [first_index-2, first_index]

    def get_stop_indices(self, strand):
        """Returns coordinates of third-to-last and last base of CDS."""
        if strand == '+':
            last_index_pair = self.indices[len(self.indices)-1]
            last_index = last_index_pair[1]
            return [last_index-2, last_index]
        elif strand == '-':
            last_index_pair = self.indices[len(self.indices)-1]
            last_index = last_index_pair[0]
            return [last_index, last_index+2]

    def sort_attributes(self):
        """Sorts indices, keeping identifiers and phases with their corresponding index pair.

        If CDS contains 'scores' they are kept sorted, too.
        """
        sort_scores = False
        length = len(self.indices)
        if length != len(self.identifier) or length != len(self.phase):
            return
        if length == len(self.score):
            sort_scores = True
        # Build a list of lists where each entry is 
        # composed of attributes
        all_attributes = []
        for i in xrange(length):
            all_attributes.append([self.indices[i][0], self.indices[i][1], 
                self.identifier[i], self.phase[i]])
            if sort_scores:
                all_attributes[i].append(self.score[i])

        # Sort that list (by first index in each index_pair)
        all_attributes.sort()
        # Repopulate the attributes
        for i in xrange(length):
            self.indices[i] = [all_attributes[i][0], all_attributes[i][1]]
            self.identifier[i] = all_attributes[i][2]
            self.phase[i] = all_attributes[i][3]
            if sort_scores:
                self.score[i] = all_attributes[i][4]

    def extract_sequence(self, seq_object, strand):
        """Returns nucleotide sequence corresponding to CDS.

        Args:
            seq_object: the actual Sequence containing the CDS. I know, I know.
            strand: either '+' or '-'
        Returns:
            a string of nucleotides (or an empty string if strand is invalid or 
            CDS has no indices)
        """
        seq = ''
        for i in xrange(len(self.indices)):
            index_pair = self.indices[i]
            subseq = seq_object.get_subseq(index_pair[0], index_pair[1])
            if subseq:
                seq += subseq
        if strand == '-':
            seq = translate.reverse_complement(seq)
        return seq

    def to_tbl(self, has_start, has_stop):
        """Returns a string containing the .tbl-formatted entry for the CDS."""
        indices = copy.deepcopy(self.indices)
        if self.strand == '+':
            phase = self.phase[0]
        else:
            phase = self.phase[-1]
        return write_tbl_entry(indices, self.strand, has_start, has_stop, "CDS", phase) 

