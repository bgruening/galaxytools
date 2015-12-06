#!/usr/bin/env python
# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

import os
import sys
from src.fasta_reader import FastaReader
from src.gff_reader import GFFReader
from src.filter_manager import FilterManager
from src.stats_manager import StatsManager

class Controller:

    def __init__(self):
        self.seqs = []
        self.removed_features = []
        self.filter_mgr = FilterManager()
        self.stats_mgr = StatsManager()

    def execute(self, args):
        """At a minimum, write a fasta, gff and tbl to output directory. Optionally do more."""
        # Verify and read fasta file
        fastapath = args.fasta
        if not os.path.isfile(fastapath):
            sys.stderr.write("Failed to find " + fastapath + ". No genome was loaded.\n")
            sys.exit()
        sys.stderr.write("Reading fasta...\n")
        self.read_fasta(fastapath)
        sys.stderr.write("Done.\n")

        # Create output directory
        out_dir = "gag_output"
        if args.out:
            out_dir = args.out
        os.system('mkdir ' + out_dir)

        # Verify and read gff file
        # This step also writes genome.ignored.gff,
        # genome.invalid.gff and genome.comments.gff
        gffpath = args.gff
        if not os.path.isfile(gffpath):
            sys.stderr.write("Failed to find " + gffpath + ". No genome was loaded.")
            return
        sys.stderr.write("Reading gff...\n")
        self.read_gff(gffpath, out_dir)
        sys.stderr.write("Done.\n")

        # Calculate stats before genome is modified
        sys.stderr.write("Calculating stats on original genome\n")
        for seq in self.seqs:
            self.stats_mgr.update_ref(seq.stats())

        # Optional annotation step
        if args.anno:
            anno_filename = args.anno
            self.annotate_from_file(anno_filename)

        # Optional step to trim sequences, subsequences or features
        if args.trim:
            trim_filename = args.trim
            self.trim_from_file(trim_filename)

        # Optional step to create start and stop codons
        if args.fix_start_stop:
            sys.stderr.write("Creating start and stop codons...\n")
            self.fix_start_stop_codons()

        # Optional step to fix terminal Ns
        if args.fix_terminal_ns:
            sys.stderr.write("Fixing terminal Ns...\n")
            self.fix_terminal_ns()

        # Optional filtering steps
        # Remove
        if args.remove_cds_shorter_than:
            min_length = args.remove_cds_shorter_than
            sys.stderr.write("Removing CDS shorter than %s...\n" % min_length)
            self.apply_filter("cds_shorter_than", min_length, "REMOVE")
        if args.remove_cds_longer_than:
            max_length = args.remove_cds_longer_than
            sys.stderr.write("Removing CDS longer than %s...\n" % max_length)
            self.apply_filter("cds_longer_than", max_length, "REMOVE")
        if args.remove_exons_shorter_than:
            min_length = args.remove_exons_shorter_than
            sys.stderr.write("Removing exons shorter than %s...\n" % min_length)
            self.apply_filter("exon_shorter_than", min_length, "REMOVE")
        if args.remove_exons_longer_than:
            max_length = args.remove_exons_longer_than
            sys.stderr.write("Removing exons longer than %s...\n" % max_length)
            self.apply_filter("exon_longer_than", max_length, "REMOVE")
        if args.remove_introns_shorter_than:
            min_length = args.remove_introns_shorter_than
            sys.stderr.write("Removing exons shorter than %s...\n" % min_length)
            self.apply_filter("intron_shorter_than", min_length, "REMOVE")
        if args.remove_introns_longer_than:
            max_length = args.remove_introns_longer_than
            sys.stderr.write("Removing exons longer than %s...\n" % max_length)
            self.apply_filter("intron_longer_than", max_length, "REMOVE")
        if args.remove_genes_shorter_than:
            min_length = args.remove_genes_shorter_than
            sys.stderr.write("Removing genes shorter than %s...\n" % min_length)
            self.apply_filter("gene_shorter_than", min_length, "REMOVE")
        if args.remove_genes_longer_than:
            max_length = args.remove_genes_longer_than
            sys.stderr.write("Removing genes longer than %s...\n" % max_length)
            self.apply_filter("gene_longer_than", max_length, "REMOVE")
        # Flag
        if args.flag_cds_shorter_than:
            min_length = args.flag_cds_shorter_than
            sys.stderr.write("Flagging CDS shorter than %s...\n" % min_length)
            self.apply_filter("cds_shorter_than", min_length, "FLAG")
        if args.flag_cds_longer_than:
            max_length = args.flag_cds_longer_than
            sys.stderr.write("Flagging CDS longer than %s...\n" % max_length)
            self.apply_filter("cds_longer_than", max_length, "FLAG")
        if args.flag_exons_shorter_than:
            min_length = args.flag_exons_shorter_than
            sys.stderr.write("Flagging exons shorter than %s...\n" % min_length)
            self.apply_filter("exon_shorter_than", min_length, "FLAG")
        if args.flag_exons_longer_than:
            max_length = args.flag_exons_longer_than
            sys.stderr.write("Flagging exons longer than %s...\n" % max_length)
            self.apply_filter("exon_longer_than", max_length, "FLAG")
        if args.flag_introns_shorter_than:
            min_length = args.flag_introns_shorter_than
            sys.stderr.write("Flagging exons shorter than %s...\n" % min_length)
            self.apply_filter("intron_shorter_than", min_length, "FLAG")
        if args.flag_introns_longer_than:
            max_length = args.flag_introns_longer_than
            sys.stderr.write("Flagging exons longer than %s...\n" % max_length)
            self.apply_filter("intron_longer_than", max_length, "FLAG")
        if args.flag_genes_shorter_than:
            min_length = args.flag_genes_shorter_than
            sys.stderr.write("Flagging genes shorter than %s...\n" % min_length)
            self.apply_filter("gene_shorter_than", min_length, "FLAG")
        if args.flag_genes_longer_than:
            max_length = args.flag_genes_longer_than
            sys.stderr.write("Flagging genes longer than %s...\n" % max_length)
            self.apply_filter("gene_longer_than", max_length, "FLAG")

        # Write fasta, gff and tbl file to output folder
        # Open files
        fasta = open(out_dir + '/genome.fasta', 'w')
        gff = open(out_dir + '/genome.gff', 'w')
        tbl = open(out_dir + '/genome.tbl', 'w')
        proteins = open(out_dir + '/genome.proteins.fasta', 'w')
        removed = open(out_dir + '/genome.removed.gff', 'w')
        stats_file = open(out_dir + '/genome.stats', 'w')

        # Calculate stats on modified genome
        sys.stderr.write("Calculating stats on modified genome\n")
        for seq in self.seqs:
            self.stats_mgr.update_alt(seq.stats())

        # Write stats file
        sys.stderr.write("Writing stats file to " + out_dir + "/ ...\n")
        for line in self.stats_mgr.summary():
            stats_file.write(line)

        # Write fasta, gff, tbl, protein fasta
        sys.stderr.write("Writing gff, tbl and fasta to " + out_dir + "/ ...\n")
        gff.write("##gff-version 3\n")
        for seq in self.seqs:
            fasta.write(seq.to_fasta())
            gff.write(seq.to_gff())
            tbl.write(seq.to_tbl())
            proteins.write(seq.to_protein_fasta())

        # Write removed.gff
        for feature in self.removed_features:
            removed.write(feature.to_gff())

        # Close files
        gff.close()
        tbl.close()
        fasta.close()
        proteins.close()
        removed.close()
        stats_file.close()

    def add_annotations_from_list(self, anno_list):
        for seq in self.seqs:
            seq.add_annotations_from_list(anno_list)

    def trim_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing trimmed.\n")
            return
        trimlist = self.read_bed_file(open(filename, 'rb'))
        if not trimlist:
            sys.stderr.write("Failed to read .bed file; nothing trimmed.\n")
            return
        else:
            self.trim_from_list(trimlist)

    def annotate_from_file(self, filename):
        if not os.path.isfile(filename):
            sys.stderr.write("Error: " + filename + " is not a file. Nothing annotated.\n")
            return
        annos = self.read_annotation_file(open(filename, 'rb'))
        if not annos:
            sys.stderr.write("Failed to read annotations from " + filename + "; no annotations added.\n")
            return
        else:
            sys.stderr.write("Adding annotations to genome ...\n")
            self.add_annotations_from_list(annos)
            sys.stderr.write("...done\n")

    def trim_from_list(self, trimlist):
        for seq in self.seqs:
            # In the case that there are multiple regions to trim in a single
            # sequence, trim from the end so indices don't get messed up
            to_trim_this_seq = [x for x in trimlist if x[0] == seq.header]
            to_trim_this_seq = sorted(to_trim_this_seq, key=lambda entry: entry[2], reverse=True)
            for entry in to_trim_this_seq:
                removed_genes = seq.trim_region(entry[1], entry[2])
                self.removed_features.extend(removed_genes)
                sys.stderr.write("Trimmed " + entry[0] + " from ")
                sys.stderr.write(str(entry[1]) + " to " + str(entry[2]) + "\n")
            self.remove_empty_features(seq)

    def get_filter_arg(self, filter_name):
        return self.filter_mgr.get_filter_arg(filter_name)
        
    def apply_filter(self, filter_name, val, filter_mode):
        for seq in self.seqs:
            self.filter_mgr.apply_filter(filter_name, val, filter_mode, seq)
            self.remove_empty_features(seq)

    def fix_terminal_ns(self):
        for seq in self.seqs:
            seq.remove_terminal_ns()
            self.remove_empty_features(seq)

    def fix_start_stop_codons(self):
        for seq in self.seqs:
            seq.create_starts_and_stops()

## Reading in files

    def read_fasta(self, line):
        reader = FastaReader()
        self.seqs = reader.read(open(line, 'r'))

    def read_gff(self, line, prefix):
        # Takes prefix b/c reader returns comments, invalids, ignored
        # and this method writes them to output files
        # That's kind of messy
        gffreader = GFFReader()
        reader = open(line, 'rb')
        genes, comments, invalids, ignored = gffreader.read_file(reader)
        for gene in genes:
            self.add_gene(gene)
        # Write comments, invalid lines and ignored features
        with open(prefix + "/genome.comments.gff", 'w') as comments_file:
            for comment in comments:
                comments_file.write(comment)
        with open(prefix + "/genome.invalid.gff", 'w') as invalid_file:
            for invalid in invalids:
                invalid_file.write(invalid)
        with open(prefix + "/genome.ignored.gff", 'w') as ignored_file:
            for item in ignored:
                ignored_file.write(item)

    def read_bed_file(self, io_buffer):
        trimlist = []
        for line in io_buffer:
            splitline = line.strip().split('\t')
            if len(splitline) != 3:
                return []
            else:
                try:
                    entry = [splitline[0], int(splitline[1]), int(splitline[2])]
                except ValueError:
                    sys.stderr.write("Error reading .bed file. Non-integer value ")
                    sys.sdterr.write("in column 2 or 3. Here is the line:\n")
                    sys.stderr.write(line)
                    return []
                trimlist.append(entry)
        return trimlist

    def read_annotation_file(self, io_buffer):
        annos = []
        for line in io_buffer:
            splitline = line.strip().split('\t')
            if len(splitline) != 3:
                return []
            else:
                annos.append(splitline)
        return annos


## Clean up

    def remove_empty_features(self, seq):
        """Removes any empty mRNAs or genes from a seq and adds them to self.removed_features."""
        self.removed_features.extend(seq.remove_empty_mrnas())
        self.removed_features.extend(seq.remove_empty_genes())
        
    def stats(self):
        if not self.seqs:
            return self.no_genome_message
        else:
            number_of_gagflags = 0
            # TODO have stats mgr handle "number of sequences"
            first_line = "Number of sequences:   " + str(len(self.seqs)) + "\n"
            sys.stderr.write("Calculating statistics on genome...\n")
            self.stats_mgr.clear_alt()
            for seq in self.seqs:
                self.stats_mgr.update_alt(seq.stats())
                number_of_gagflags += seq.number_of_gagflags()
            last_line = "(" + str(number_of_gagflags) + " features flagged)\n"
            return first_line + self.stats_mgr.summary() + last_line

## Utility methods

    def add_gene(self, gene):
        for seq in self.seqs:
            if seq.header == gene.seq_name:
                seq.add_gene(gene)

    def get_locus_tag(self):
        locus_tag = ""
        for seq in self.seqs:
            if locus_tag:
                break
            else:
                locus_tag = seq.get_locus_tag()
        return locus_tag
    
    def remove_from_list(self, bad_list):
        # First remove any seqs on the list
        to_remove = []
        for seq in self.seqs:
            if seq.header in bad_list:
                to_remove.append(seq)
        if to_remove:
            for seq in to_remove:
                self.seqs.remove(seq)
                sys.stderr.write("Warning: removing seq " + seq.header + ".\n")
                sys.stderr.write("You must reload genome to get this sequence back.\n")
            self.removed_features.extend(to_remove)
        # Now pass the list down to each seq
        for seq in self.seqs:
            removed_from_seq = seq.remove_from_list(bad_list)
            self.removed_features.extend(removed_from_seq)

    def contains_mrna(self, mrna_id):
        for seq in self.seqs:
            if seq.contains_mrna(mrna_id):
                return True
        return False

    def contains_gene(self, gene_id):
        for seq in self.seqs:
            if seq.contains_gene(gene_id):
                return True
        return False
