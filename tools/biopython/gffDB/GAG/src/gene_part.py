#!/usr/bin/env python

import math
import src.translator as translate

class GenePart:
    def __init__(self, feature_type=None, identifier=None,\
                 indices=None, score=None, strand='+', parent_id=None):
        self.feature_type = feature_type
        self.identifier = []
        if identifier is not None:
            self.identifier.append(identifier)
        self.indices = []
        if indices is not None:
            self.indices.append(indices)
        self.score = []
        if score is not None:
            self.score.append(score)
        self.strand = strand # Defauts to positive strand (?)
        self.parent_id = parent_id
        self.annotations = []

    def __str__(self):
        """Returns string representation of a GenePart.

        String contains the type of feature and the identifier for its first segment.
        """
        result = self.feature_type + " (first ID="
        result += str(self.identifier[0])+ ")"
        return result

    def add_indices(self, ind):
        """Adds an index pair to the GenePart's indices list.

        Also sorts the indices afterward in case a pair is added out of order.
        """
        if isinstance(ind, list) and len(ind) is 2:
            self.indices.append(ind)
        else:
            raise ValueError()

    def add_identifier(self, identifier):
        """Adds an identifier to the GenePart's list of identifiers."""
        self.identifier.append(identifier)

    def add_score(self, score):
        """Adds a score to the GenePart's list of scores."""
        self.score.append(score)

    def sort_attributes(self):
        """Sorts indices, keeping identifiers with their corresponding index pair.

        If GenePart contains 'scores' they are kept sorted, too.
        """
        sort_scores = False
        length = len(self.indices)
        if length != len(self.identifier):
            return
        if length == len(self.score):
            sort_scores = True
        # Build a list of lists where each entry is 
        # composed of attributes
        all_attributes = []
        for i in xrange(length):
            all_attributes.append([self.indices[i][0], self.indices[i][1], 
                self.identifier[i]])
            if sort_scores:
                all_attributes[i].append(self.score[i])

        # Sort that list (by first index in each index_pair)
        all_attributes.sort()
        # Repopulate the attributes
        for i in xrange(length):
            self.indices[i] = [all_attributes[i][0], all_attributes[i][1]]
            self.identifier[i] = all_attributes[i][2]
            if sort_scores:
                self.score[i] = all_attributes[i][3]

    def add_annotation(self, key, value):
        """Adds an annotation key, value pair to the GenePart.

        Args:
            key: a string indicating the type of annotation (e.g. Dbxref, gagflag)
            value: a string representing the content of the annotation
        """
        self.annotations.append([key, value])
    
    def gagflagged(self):
        """Returns a boolean indicating whether the GenePart contains a 'gagflag' annotation."""
        for anno in self.annotations:
            if anno[0] == "gag_flag":
                return True
        return False

    def length(self):
        """Returns the entire length of a feature."""
        length = 0
        for index_pair in self.indices:
            length += length_of_segment(index_pair)
        return length

    # used by .to_gff
    def get_score(self, i):
        """Returns the 'score' entry corresponding to a given index pair of a feature.

        Returns '.' if no score entry is present; this is the default value in .gff format.
        """
        if self.score and len(self.score) > i:
            return self.score[i]
        else:
            return "."

    # used by .to_gff
    # (CDS overrides this method)
    def get_phase(self, i):
        """Returns the phase corresponding to a given index pair of a feature.

        This method exists to be overridden by CDS; it simply returns '.'
        for other features as this is the default value in a .gff file.
        """
        return "."

    def adjust_indices(self, n, start_index=1):
        """Increments indices of GenePart.

        Optionally, only indices occurring after start_index are incremented.

        Args:
            n: integer by which to increment indices
            start_index: optional coordinate before which no indices will be changed.
        """
        for i, index_pair in enumerate(self.indices):
            if index_pair[0] >= start_index:
                self.indices[i] = adjust_index_pair(self.indices[i], n)
            elif index_pair[1] >= start_index:
                self.indices[i][1] += n

    def generate_attribute_entry(self, i):
        """Returns a string representing a GenePart's .gff attribute entry.

        Includes mandatory 'ID', 'Parent' fields as well as annotations.
        Args:
            i: index of the row of .gff data to be returned.
        """
        if len(self.identifier) <= i or self.parent_id is None:
            return None
        entry = "ID=" + str(self.identifier[i]) + ";"
        entry += "Parent=" + str(self.parent_id)
        for annot in self.annotations:
            entry += ';'+annot[0]+'='+annot[1]
        entry += '\n'
        return entry

    def to_gff(self, seq_name, source):
        """Returns a string containing the .gff representation of a GenePart."""
        result = ""
        for i in range(len(self.indices)):
            result += seq_name + "\t" + source + "\t"
            result += self.feature_type + "\t" + str(self.indices[i][0])
            result += "\t" + str(self.indices[i][1]) + "\t"
            result += str(self.get_score(i)) + "\t" + self.strand + "\t"
            result += str(self.get_phase(i)) + "\t"
            result += self.generate_attribute_entry(i)
        return result

######################### Utility Functions ##########################################

def length_of_segment(index_pair):
    """Returns the length in base pairs between two indices (inclusive)."""
    return math.fabs(index_pair[1] - index_pair[0]) + 1

def adjust_index_pair(index_pair, increment):
    """Returns pair of indices incremented by given number."""
    return [i + increment for i in index_pair]

def get_reversed_indices(indices):
    """Returns reversed list of indices.

    Each pair in the list is reversed, and the order of the
    pairs is also reversed.
    """
    indices.reverse()
    [ind.reverse() for ind in indices]
    return indices

def one_line_indices_entry(indices, has_start, has_stop, feature_type):
    """Returns .tbl formatted entry for features with only one pair of coordinates."""
    output = ""
    if not has_start:
        output += "<"
    output += str(indices[0]) + "\t"
    if not has_stop:
        output += ">"
    output += str(indices[1]) + "\t"
    output += feature_type+"\n"
    return output

def write_tbl_entry(indices, strand, has_start, has_stop, feature_type, phase=0):
    """Returns .tbl-formatted string for a GenePart.

    Args:
        indices: a list of lists, each holding an index pair
        strand: either '+' or '-'
        has_start: a boolean indicating whether the feature has a start codon
        has_stop: a boolean indicating whether the feature has a stop codon
        feature_type: tbl feature type (i.e. CDS, mRNA, tRNA, etc)
        phase: optional argument indicating the phase of a CDS feature if not 0
    """
    output = ""
    if strand == "-":
        indices = get_reversed_indices(indices)
    if len(indices) == 1:
        output += one_line_indices_entry(indices[0], has_start, has_stop, feature_type)
    else:
        # Write first line of coordinates
        if not has_start:
            output += "<"
        output += str(indices[0][0]) + "\t" + str(indices[0][1]) + "\t" 
        output += feature_type+"\n"
        # Write middle lines
        for index_pair in indices[1:-1]:
            output += str(index_pair[0]) + "\t" + str(index_pair[1]) + "\n"
        # Write last line of coordinates
        output += str(indices[-1][0]) + "\t"
        if not has_stop:
            output += ">"
        output += str(indices[-1][1]) + "\n"
    if feature_type == "CDS":
        output += "\t\t\tcodon_start\t" + str(phase+1) + "\n"
    return output

