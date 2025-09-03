#!/usr/bin/env python

import math
import sys

def length_of_segment(index_pair):
    return math.fabs(index_pair[1] - index_pair[0]) + 1

class Gene:

    def __init__(self, seq_name, source, indices, strand, identifier, name="", annotations=None, score=None):
        self.seq_name = seq_name
        self.source = source
        self.indices = indices
        self.score = score
        self.strand = strand
        self.identifier = identifier
        self.name = name
        self.mrnas = []
        self.removed_mrnas = []
        if not annotations:
            self.annotations = {}
        else:
            self.annotations = annotations
        self.death_flagged = False

    def __str__(self):
        """Returns string representation of a gene.

        String contains the gene's identifier, the sequence it is on,
        and the number of mRNAs it contains.
        """
        result = "Gene (ID=" + str(self.identifier) 
        result += ", seq_name=" + self.seq_name
        result += ") containing " + str(len(self.mrnas))
        result += " mrnas"
        return result

    def get_mrna_ids(self):
        """Returns a list of the identifiers of all mRNAs contained on the gene."""
        result = []
        for mrna in self.mrnas:
            result.append(mrna.identifier)
        return result
    
    def remove_mrna(self, mrna_id):
        to_remove = None
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id:
                to_remove = mrna
        if to_remove:
            self.mrnas.remove(to_remove)
            self.removed_mrnas.append(to_remove)
            return True
        return False # Return false if mrna wasn't removed
        
    def remove_mrnas_from_list(self, bad_mrnas):
        to_remove = []
        for mrna in self.mrnas:
            if mrna.identifier in bad_mrnas:
                to_remove.append(mrna)
        if to_remove:
            for mrna in to_remove:
                self.mrnas.remove(mrna)
                sys.stderr.write("Removed mrna " + mrna.identifier + "\n")
            self.removed_mrnas.extend(to_remove)
        return to_remove
    
    def remove_empty_mrnas(self):
        """Removes mRNAs with no exon or CDS.
        """
        to_remove = []
        for mrna in self.mrnas:
            if not mrna.cds and not mrna.exon:
                to_remove.append(mrna)
                sys.stderr.write("Removed empty mrna " + mrna.identifier + "\n")
            elif mrna.rna_type == "mRNA" and not mrna.cds:
                to_remove.append(mrna)
                sys.stderr.write("Removed mrna " + mrna.identifier + " with no cds\n")
            elif not mrna.exon:
                to_remove.append(mrna)
                sys.stderr.write("Removed mrna " + mrna.identifier + " with no exon\n")
        if to_remove:
            for mrna in to_remove:
                self.mrnas.remove(mrna)
                
            self.removed_mrnas.extend(to_remove)
        return to_remove

    def add_annotation(self, key, value):
        """Adds an annotation key, value pair to the gene.

        Args:
            key: a string indicating the type of annotation (e.g. Dbxref, gagflag)
            value: a string representing the content of the annotation
        """
        if key in self.annotations.keys():
            self.annotations[key].append(value)
        else:
            self.annotations[key] = [value]

    def add_mrna_annotation(self, mrna_id, key, value):
        """Adds annotation key, value pair to specified mrna.

        Does nothing if mrna not found.
        Args:
            mrna_id: the identifier of the mrna to which the annotation belongs
            key: a string indicating the type of annotation (e.g. Dbxref, gagflag)
            value: a string representing the content of the annotation
        """
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id:
                mrna.add_annotation(key, value)
        
    def length(self):
        """Returns the length of the gene."""
        return length_of_segment(self.indices)

    def gagflagged(self):
        """Returns a boolean indicating whether the gene itself contains a 'gagflag' annotation."""
        return 'gag_flag' in self.annotations.keys()

    def number_of_gagflags(self):
        """Returns the number of gagflagged features on the gene, including the gene itself."""
        total = 0
        for mrna in self.mrnas:
            total += mrna.number_of_gagflags()
        if self.gagflagged():
            total += 1
        return total

    def get_score(self):
        """Returns the gene's score (column 6 in .gff files).

        If no score is present, returns '.' (the default .gff value).
        """
        if self.score:
            return self.score
        else:
            return '.'

    def get_longest_exon(self):
        """Returns length of longest exon contained on gene."""
        longest = 0
        for mrna in self.mrnas:
            length = mrna.get_longest_exon()
            if length > longest:
                longest = length
        return longest

    def get_shortest_exon(self):
        """Returns length of shortest exon contained on gene."""
        shortest = None
        for mrna in self.mrnas:
            length = mrna.get_shortest_exon()
            if length == 0:
                continue
            if shortest == None or length < shortest:
                shortest = length
        if shortest == None:
            return 0
        return shortest

    def get_total_exon_length(self):
        """Returns sum of all child exon lengths."""
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_total_exon_length()
        return total
    
    def get_num_exons(self):
        """Returns number of exons contained on gene."""
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_num_exons()
        return total

    def get_longest_intron(self):
        """Returns length of longest intron contained on gene."""
        longest = 0
        for mrna in self.mrnas:
            length = mrna.get_longest_intron()
            if length > longest:
                longest = length
        return longest

    def get_shortest_intron(self):
        """Returns length of shortest intron contained on gene."""
        shortest = None
        for mrna in self.mrnas:
            length = mrna.get_shortest_intron()
            if length == 0:
                continue
            if shortest == None or length < shortest:
                shortest = length
        if shortest == None:
            return 0
        return shortest

    def get_total_intron_length(self):
        """Returns sum of all child intron lengths."""
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_total_intron_length()
        return total

    def get_num_introns(self):
        """Returns number of introns contained on gene."""
        total = 0
        for mrna in self.mrnas:
            total += mrna.get_num_introns()
        return total
    
    def create_starts_and_stops(self, seq_object):
        """Creates start and stop codons on child mRNAs.

        Args:
            seq_object: the actual Sequence containing the gene. I know, I know.
        """
        for mrna in self.mrnas:
            mrna.create_start_and_stop_if_necessary(seq_object, self.strand)

    def adjust_indices(self, n, start_index=1):
        """Adds 'n' to both indices, checking to ensure that they fall after an optional start index"""
        if self.indices[0] >= start_index:
            self.indices = [i + n for i in self.indices]
        elif self.indices[1] >= start_index:
            self.indices[1] += n
        for mrna in self.mrnas:
            mrna.adjust_indices(n, start_index)

    def get_partial_info(self):
        """Returns a dictionary containing counts for complete/incomplete CDSs."""
        results = {"complete": 0, "start_no_stop": 0, "stop_no_start": 0, "no_stop_no_start": 0}
        for mrna in self.mrnas:
            if mrna.has_start():
                if mrna.has_stop():
                    results["complete"] += 1
                else:
                    results["start_no_stop"] += 1
            else:
                # No start ...
                if mrna.has_stop():
                    results["stop_no_start"] += 1
                else:
                    results["no_stop_no_start"] += 1
        return results

    def remove_mrnas_with_internal_stops(self, seq_helper):
        """Removes child mRNAs that contain internal stop codons."""
        for mrna in self.mrnas:
            if seq_helper.mrna_contains_internal_stop(mrna):
                mrna.death_flagged = True
        self.mrnas = [m for m in self.mrnas if not m.death_flagged]

    def contains_mrna(self, mrna_id):
        """Returns a boolean indicating whether gene contains an mRNA with the given id."""
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id:
                return True
        return False

    def cds_to_gff(self, seq_id, mrna_id):
        """Returns a string containing the gff representation of a CDS.

        If the CDS is not contained on the gene, returns an empty string.

        Args:
            seq_id: the 'header' field of the current seq
            mrna_id: the id of the mRNA containing the CDS to be written.
        """
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id and mrna.cds:
                return mrna.cds_to_gff(seq_id, self.source)
        return ""

    def cds_to_tbl(self, mrna_id):
        """Returns a string containing the tbl representation of a CDS.

        If the CDS is not contained on the gene, returns an empty string.

        Args:
            mrna_id: the id of the mRNA containing the CDS to be written.
        """
        for mrna in self.mrnas:
            if mrna.identifier == mrna_id and mrna.cds:
                return mrna.cds_to_tbl()
        return ""

    def to_mrna_fasta(self, seq_helper):
        """Returns a string containing fasta entries for the gene's mRNAs."""
        result = ""
        for mrna in self.mrnas:
            result += seq_helper.mrna_to_fasta(mrna)
        return result

    def to_cds_fasta(self, seq_helper):
        """Returns a string containing fasta entries for the gene's CDSs."""
        result = ""
        for mrna in self.mrnas:
            result += seq_helper.mrna_to_cds_fasta(mrna)
        return result

    def to_protein_fasta(self, seq_helper):
        """Returns a string containing fasta entries for the gene's proteins."""
        result = ""
        for mrna in self.mrnas:
            result += seq_helper.mrna_to_protein_fasta(mrna)
        return result

    def to_gff(self, removed_features=False):
        """Returns a string in .gff format of the gene and its child features."""
        result = self.seq_name + "\t" + self.source + "\t"
        result += 'gene' + "\t" + str(self.indices[0]) + "\t"
        result += str(self.indices[1]) + "\t" + self.get_score()
        result += "\t" + self.strand + "\t" + "." + "\t"
        result += "ID=" + str(self.identifier)
        if self.name:
            result += ";Name=" + self.name
        for key in self.annotations.keys():
            result += ';' + key + "="
            result += ','.join(self.annotations[key])
        result += '\n'
        for mrna in self.mrnas:
            result += mrna.to_gff()
        # Now write the removed features if they want them
        if removed_features:
            for mrna in self.removed_mrnas:
                result += mrna.to_gff()
        return result
    
    # Outputs only removed mrnas
    def removed_to_gff(self):
        """Returns a string in .gff format of the gene and its child features."""
        result = ""
        for mrna in self.removed_mrnas:
            result += mrna.to_gff()
        return result

    def to_tbl(self):
        """Returns a string in .tbl format of the gene and its child features."""
        if self.strand == "-":
            indices = [self.indices[1], self.indices[0]]
        else:
            indices = self.indices
        # Check if there's an mRNA with a no start/stop
        has_start = True
        has_stop = True
        for mrna in self.mrnas:
            if not mrna.has_start():
                has_start = False
            if not mrna.has_stop():
                has_stop = False
        output = ""
        if not has_start:
            output += "<"
        output += str(indices[0]) + "\t"
        if not has_stop:
            output += ">"
        output += str(indices[1]) + "\t" + "gene\n"
        if self.name:
            output += "\t\t\tgene\t" + self.name + "\n"
        output += "\t\t\tlocus_tag\t" + self.identifier + "\n"
        for mrna in self.mrnas:
            output += mrna.to_tbl()
        return output



