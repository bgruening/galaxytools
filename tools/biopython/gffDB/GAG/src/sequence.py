#!/usr/bin/env python

import sys
from src.seq_helper import SeqHelper

class Sequence:

    def __init__(self, header="", bases=""):
        self.header = header
        self.bases = bases
        self.genes = []
        self.removed_genes = []

    def __str__(self):
        result = "Sequence " + self.header
        result += " of length " + str(len(self.bases))
        result += " containing "
        result += str(len(self.genes))
        result += " genes\n"
        return result

    def add_gene(self, gene):
        self.genes.append(gene)

    def contains_gene(self, gene_id):
        for gene in self.genes:
            if gene.identifier == gene_id:
                return True
        return False

    def contains_mrna(self, mrna_id):
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.identifier == mrna_id:
                    return True
        return False
    
    def remove_gene(self, gene_id):
        to_remove = None
        for gene in self.genes:
            if gene.identifier == gene_id:
                to_remove = gene
        if to_remove:
            self.genes.remove(to_remove)
            self.removed_genes.append(to_remove)
            return True
        return False # Return false if gene wasn't removed

    def remove_from_list(self, bad_list):
        removed_features = []
        removed_genes = self.remove_genes_from_list(bad_list)
        removed_mrnas = self.remove_mrnas_from_list(bad_list)
        removed_features.extend(removed_genes)
        removed_features.extend(removed_mrnas)
        return removed_mrnas

    def remove_genes_from_list(self, bad_genes):
        to_remove = []
        for gene in self.genes:
            if gene.identifier in bad_genes:
                to_remove.append(gene)
        if to_remove:
            for gene in to_remove:
                self.genes.remove(gene)
                sys.stderr.write("Removed gene " + gene.identifier + "\n")
            self.removed_genes.extend(to_remove)
        return to_remove

    def remove_mrnas_from_list(self, bad_mrnas):
        removed_mrnas = []
        for gene in self.genes:
            removed_from_gene = gene.remove_mrnas_from_list(bad_mrnas)
            removed_mrnas.extend(removed_from_gene)
        return removed_mrnas
    
    def remove_empty_genes(self):
        """Removes any gene containing no mRNAs; returns a list of removed genes."""
        to_remove = []
        for gene in self.genes:
            if not gene.mrnas:
                to_remove.append(gene)
        if to_remove:
            for gene in to_remove:
                self.genes.remove(gene)
                sys.stderr.write("Removed empty gene " + gene.identifier + "\n")
            self.removed_genes.extend(to_remove)
        return to_remove
    
    def remove_empty_mrnas(self):
        removed_mrnas = []
        for gene in self.genes:
            removed_from_gene = gene.remove_empty_mrnas()
            removed_mrnas.extend(removed_from_gene)
        return removed_mrnas

    def add_annotations_from_list(self, anno_list):
        for gene in self.genes:
            for anno in anno_list:
                if gene.identifier == anno[0] and anno[1] == "name":
                    gene.name = anno[2]
                if gene.contains_mrna(anno[0]):
                    gene.add_mrna_annotation(anno[0], anno[1], anno[2])

    def number_of_gagflags(self):
        total = 0
        for gene in self.genes:
            total += gene.number_of_gagflags()
        return total
        
    def get_gene_ids(self):
        result = []
        for gene in self.genes:
            result.append(gene.identifier)
        return result

    def get_mrna_ids(self):
        result = []
        for gene in self.genes:
            result.extend(gene.get_mrna_ids())
        return result

    def get_locus_tag(self):
        for gene in self.genes:
            gene_id = str(gene.identifier)
            locus_tag = gene_id.split('_')[0]
            return locus_tag
        return ""

    def to_fasta(self):
        result = '>' + self.header + '\n'
        result += self.bases + '\n'
        return result

    def remove_terminal_ns(self):
        # Remove any Ns at the beginning of the sequence
        initial_ns = self.how_many_Ns_forward(1)
        if initial_ns:
            self.trim_region(1, initial_ns)
        # Remove any Ns at the end of the sequence
        length = len(self.bases)
        terminal_ns = self.how_many_Ns_backward(length)
        if terminal_ns:
            self.trim_region(length-terminal_ns+1, length)

    # Given a position in the sequence, returns the number of Ns 
    # from that position forward 
    # (returns 0 if the base at that position is not N)
    def how_many_Ns_forward(self, position):
        index = position-1
        if self.bases[index] != 'N' and self.bases[index] != 'n':
            return 0
        else:
            count = 1
            index += 1
            for base in self.bases[index:]:
                if base != 'N' and base != 'n':
                    return count
                else:
                    count += 1
            return count

    # Given a position in the fasta, returns the number of Ns 
    # from that position backward 
    # (returns 0 if the base at that position is not N)
    def how_many_Ns_backward(self, position):
        index = position-1
        if self.bases[index] != 'N' and self.bases[index] != 'n':
            return 0
        else:
            count = 1
            index -= 1
            for base in self.bases[index::-1]:
                if base != 'N' and base != 'n':
                    return count
                else:
                    count += 1
            return count

    def trim_region(self, start, stop):
        """Remove bases from start to stop; remove and return affected genes"""
        if stop > len(self.bases):
            sys.stderr.write("Sequence.trim called on sequence that is too short;"+\
                    " doing nothing.\n")
            return
        # Remove any genes that are overlap the trimmed region
        genes_to_remove = [g for g in self.genes if overlap([start, stop], g.indices)]
        self.genes = [g for g in self.genes if g not in genes_to_remove]
        # Remove bases from sequence
        self.bases = self.bases[:start-1] + self.bases[stop:]
        # Adjust indices of remaining genes
        bases_removed = stop - start + 1
        [g.adjust_indices(-bases_removed, start) for g in self.genes]
        return genes_to_remove
        
    def get_subseq(self, start=1, stop=None):
        if not stop:
            stop = len(self.bases)
        if stop > len(self.bases):
            return ""
        return self.bases[start-1:stop]

    def remove_mrnas_with_internal_stops(self):
        helper = SeqHelper(self.bases)
        for gene in self.genes:
            gene.remove_mrnas_with_internal_stops(helper)
            if not gene.mrnas:
                gene.death_flagged = True
        self.genes = [g for g in self.genes if not g.death_flagged]

    def create_starts_and_stops(self):
        for gene in self.genes:
            gene.create_starts_and_stops(self)
    
    def get_contained_genes(self):
        contained = []
        for i, a in enumerate(self.genes):
            for b in self.genes[i+1:]:
                # Skip duplicate indices
                if a.indices == b.indices:
                    continue
                # Check if a contains b
                if a.indices[0] <= b.indices[0] and a.indices[1] >= b.indices[1] and not b in contained:
                    contained.append(b)
                # Check if b contains a
                elif b.indices[0] <= a.indices[0] and b.indices[1] >= a.indices[1] and not a in contained:
                    contained.append(a)
        return contained
    
    def get_overlapping_genes(self):
        contained = []
        for i, a in enumerate(self.genes):
            for b in self.genes[i+1:]:
                if a.indices[1] >= b.indices[0] and a.indices[0] <= b.indices[1]:
                    if not a in contained:
                        contained.append(a)
                    if not b in contained:
                        contained.append(b)
        return contained

    def cds_to_gff(self, mrna_id):
        for gene in self.genes:
            if gene.contains_mrna(mrna_id):
                return gene.cds_to_gff(self.header, mrna_id)
        return "CDS not found."

    def cds_to_tbl(self, mrna_id):
        for gene in self.genes:
            if gene.contains_mrna(mrna_id):
                return gene.cds_to_tbl(mrna_id)
        return "CDS not found."

    def to_tbl(self):
        result = ">Feature " + self.header + "\n"
        result += "1\t" + str(len(self.bases)) + "\tREFERENCE\n"
        result += "\t\t\tPBARC\t12345\n"
        for gene in self.genes:
            result += gene.to_tbl()
        return result

    def to_mrna_fasta(self):
        helper = SeqHelper(self.bases)
        result = ""
        for gene in self.genes:
            result += gene.to_mrna_fasta(helper)
        return result

    def to_cds_fasta(self):
        helper = SeqHelper(self.bases)
        result = ""
        for gene in self.genes:
            result += gene.to_cds_fasta(helper)
        return result

    def to_protein_fasta(self):
        helper = SeqHelper(self.bases)
        result = ""
        for gene in self.genes:
            result += gene.to_protein_fasta(helper)
        return result

    def to_gff(self):
        result = ""
        for gene in self.genes:
            result += gene.to_gff()
        return result
    
    def removed_to_gff(self):
        result = ""
        # Write alive genes' removed mrnas
        for gene in self.genes:
            result += gene.removed_to_gff()
        # Write all dead genes' mrnas
        for gene in self.removed_genes:
            result += gene.to_gff(True)
        return result

###################################################################################################
# Statsy type stuff

    def get_num_mrna(self):
        count = 0
        for gene in self.genes:
            count += len(gene.mrnas)
        return count

    def get_num_exons(self):
        count = 0
        for gene in self.genes:
            count += gene.get_num_exons()
        return count
        
    def get_num_cds(self):
        count = 0
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.cds:
                    count += 1
        return count

    def get_cds_partial_info(self):
        results = {"CDS: complete": 0, "CDS: start, no stop": 0,\
                "CDS: stop, no start": 0, "CDS: no stop, no start": 0}
        for gene in self.genes:
            partial_info = gene.get_partial_info()
            results["CDS: complete"] += partial_info["complete"]
            results["CDS: start, no stop"] += partial_info["start_no_stop"]
            results["CDS: stop, no start"] += partial_info["stop_no_start"]
            results["CDS: no stop, no start"] += partial_info["no_stop_no_start"]
        return results
        
    def get_longest_gene(self):
        length = 0
        for gene in self.genes:
            if gene.length() > length:
                length = gene.length()
        return length
        
    def get_longest_mrna(self):
        length = 0
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.length() > length:
                    length = mrna.length()
        return length

    def get_longest_exon(self):
        longest = 0
        for gene in self.genes:
            length = gene.get_longest_exon()
            if length > longest:
                longest = length
        return longest

    def get_shortest_exon(self):
        shortest = None
        for gene in self.genes:
            length = gene.get_shortest_exon()
            if length == 0:
                continue
            if shortest == None or length < shortest:
                shortest = length
        if shortest == None:
            return 0
        return shortest

    def get_total_exon_length(self):
        total = 0
        for gene in self.genes:
            total += gene.get_total_exon_length()
        return total
        
    def get_longest_intron(self):
        longest = 0
        for gene in self.genes:
            length = gene.get_longest_intron()
            if length > longest:
                longest = length
        return longest

    def get_shortest_intron(self):
        shortest = None
        for gene in self.genes:
            length = gene.get_shortest_intron()
            if length == 0:
                continue
            if shortest == None or length < shortest:
                shortest = length
        if shortest == None:
            return 0
        return shortest

    def get_total_intron_length(self):
        total = 0
        for gene in self.genes:
            total += gene.get_total_intron_length()
        return total

    def get_num_introns(self):
        total = 0
        for gene in self.genes:
            total += gene.get_num_introns()
        return total
        
    def get_longest_cds(self):
        length = 0
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.cds and mrna.cds.length() > length:
                    length = mrna.cds.length()
        return length
        
    def get_shortest_gene(self):
        length = 0
        shortest = None
        for gene in self.genes:
            if gene.length() < length or shortest == None:
                length = gene.length()
                shortest = gene
        return length
        
    def get_shortest_mrna(self):
        length = 0
        shortest = None
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.length() < length or shortest == None:
                    length = mrna.length()
                    shortest = mrna
        return length
        
    def get_shortest_cds(self):
        length = 0
        shortest = None
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.cds and (mrna.cds.length() < length or shortest == None):
                    length = mrna.cds.length()
                    shortest = mrna.cds
        return length

    def get_total_gene_length(self):
        length = 0
        for gene in self.genes:
            length += gene.length()
        return length

    def get_total_mrna_length(self):
        length = 0
        for gene in self.genes:
            for mrna in gene.mrnas:
                length += mrna.length()
        return length
        
    def get_total_cds_length(self):
        length = 0
        for gene in self.genes:
            for mrna in gene.mrnas:
                if mrna.cds:
                    length += mrna.cds.length()
        return length
        
    def stats(self):
        stats = dict()
        cds_partial_info = self.get_cds_partial_info()
        
        stats["Total sequence length"] = len(self.bases)
        stats["Number of genes"] = len(self.genes)
        stats["Number of mRNAs"] = int(self.get_num_mrna())
        stats["Number of exons"] = int(self.get_num_exons())
        stats["Number of introns"] = int(self.get_num_introns())
        stats["Number of CDS"] = int(self.get_num_cds())
        stats["Overlapping genes"] = len(self.get_overlapping_genes())
        stats["Contained genes"] = len(self.get_contained_genes())
        stats["CDS: complete"] = int(cds_partial_info["CDS: complete"])
        stats["CDS: start, no stop"] = int(cds_partial_info["CDS: start, no stop"])
        stats["CDS: stop, no start"] = int(cds_partial_info["CDS: stop, no start"])
        stats["CDS: no stop, no start"] = int(cds_partial_info["CDS: no stop, no start"])
        stats["Longest gene"] = int(self.get_longest_gene())
        stats["Longest mRNA"] = int(self.get_longest_mrna())
        stats["Longest exon"] = int(self.get_longest_exon())
        stats["Longest intron"] = int(self.get_longest_intron())
        stats["Longest CDS"] = int(self.get_longest_cds())
        stats["Shortest gene"] = int(self.get_shortest_gene())
        stats["Shortest mRNA"] = int(self.get_shortest_mrna())
        stats["Shortest exon"] = int(self.get_shortest_exon())
        stats["Shortest intron"] = int(self.get_shortest_intron())
        stats["Shortest CDS"] = int(self.get_shortest_cds())
        stats["Total gene length"] = int(self.get_total_gene_length())
        stats["Total mRNA length"] = int(self.get_total_mrna_length())
        stats["Total exon length"] = int(self.get_total_exon_length())
        stats["Total intron length"] = int(self.get_total_intron_length())
        stats["Total CDS length"] = int(self.get_total_cds_length())
        
        return stats

def overlap(indices1, indices2):
    """Returns a boolean indicating whether two pairs of indices overlap."""
    if not (len(indices1) == 2 and len(indices2) ==2):
        return False
    if indices1[0] >= indices2[0] and indices1[0] <= indices2[1]:
        return True
    elif indices1[1] >= indices2[0] and indices1[1] <= indices2[1]:
        return True
    else:
        return False
