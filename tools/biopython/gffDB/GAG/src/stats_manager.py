#!/usr/bin/env python


### TODO TODO TODO
### ok, update should take a list of seqs instead of a single seq
### then it can calculate stuff like longest seq, number of seqs, etc.

class StatsManager:

    increment_stats = ["Total sequence length", "Number of genes", "Number of mRNAs", "Number of exons", "Number of introns", "Number of CDS",\
            "Overlapping genes", "Contained genes", "CDS: complete", "CDS: start, no stop", "CDS: stop, no start", "CDS: no stop, no start",\
            "Total gene length", "Total mRNA length", "Total exon length",\
            "Total intron length", "Total CDS length"]    
    min_stats = ["Shortest gene", "Shortest mRNA", "Shortest exon", "Shortest intron", "Shortest CDS"]
    max_stats = ["Longest gene", "Longest mRNA", "Longest exon", "Longest intron", "Longest CDS"]
    calc_stats_formulae = {"mean gene length": ["Total gene length", "Number of genes"],\
            "mean mRNA length": ["Total mRNA length", "Number of mRNAs"],\
            "mean exon length": ["Total exon length", "Number of exons"],\
            "mean intron length": ["Total intron length", "Number of introns"],\
            "mean CDS length": ["Total CDS length", "Number of CDS"],\
            "% of genome covered by genes": ["Total gene length", "Total sequence length"],\
            "% of genome covered by CDS": ["Total CDS length", "Total sequence length"],\
            "mean mRNAs per gene": ["Number of mRNAs", "Number of genes"], "mean exons per mRNA": ["Number of exons", "Number of mRNAs"],\
            "mean introns per mRNA": ["Number of introns", "Number of mRNAs"]}
    calc_stats = ["mean gene length", "mean mRNA length", "mean exon length", "mean intron length",\
            "mean CDS length", "% of genome covered by genes", "% of genome covered by CDS",\
            "mean mRNAs per gene", "mean exons per mRNA", "mean introns per mRNA"]

    def __init__(self):
        self.ref_stats = {}
        self.alt_stats = {}
        self.initialize_dict(self.ref_stats)
        self.initialize_dict(self.alt_stats)

    def initialize_dict(self, d):
        for stat in self.increment_stats + self.min_stats + self.max_stats + self.calc_stats:
            d[stat] = 0

    def alt_is_empty(self):
        for key in self.ref_stats.keys():
            if self.ref_stats[key] != self.alt_stats[key]:
                return False
        return True

    def clear_alt(self):
        self.initialize_dict(self.alt_stats)

    def clear_all(self):
        self.initialize_dict(self.ref_stats)
        self.initialize_dict(self.alt_stats)

    def update_ref(self, stats):
        self.update_stats(self.ref_stats, stats)

    def update_alt(self, stats):
        self.update_stats(self.alt_stats, stats)

    def update_stats(self, old, new):
        if not validate_dicts(old, new):
            return
        for stat in self.increment_stats:
            old[stat] += new[stat]
        for stat in self.min_stats:
            if old[stat] == 0:
                old[stat] = new[stat]
            elif new[stat] < old[stat] and new[stat] != 0:
                old[stat] = new[stat]
        for stat in self.max_stats:
            if new[stat] > old[stat]:
                old[stat] = new[stat]

    def calculate_stat(self, stat):
        dividend_key = self.calc_stats_formulae[stat][0]
        divisor_key = self.calc_stats_formulae[stat][1]
        # Calculate for reference genome
        dividend = self.ref_stats[dividend_key]
        divisor = self.ref_stats[divisor_key]
        if divisor == 0:
            self.ref_stats[stat] = 0
        else:
            raw_value = float(dividend) / float(divisor)
            if "%" in stat:
                self.ref_stats[stat] = format_percent(raw_value)
            else:
                self.ref_stats[stat] = int(round(raw_value))
        # Calculate for modified genome
        dividend = self.alt_stats[dividend_key]
        divisor = self.alt_stats[divisor_key]
        if divisor == 0:
            self.alt_stats[stat] = 0
        else:
            raw_value = float(dividend) / float(divisor)
            if "%" in stat:
                self.alt_stats[stat] = format_percent(raw_value)
            else:
                self.alt_stats[stat] = int(round(raw_value))

    def summary(self):
        for stat in self.calc_stats:
            self.calculate_stat(stat)
        stats_order = [key for keys in [self.increment_stats, self.min_stats, self.max_stats, self.calc_stats] for key in keys]
        if self.alt_is_empty():
            return format_columns(["Genome"], stats_order, [self.ref_stats], 5)
        else:
            return format_columns(["Reference Genome", "Modified Genome"], stats_order, [self.ref_stats, self.alt_stats], 5)

## UTILITY FUNCTIONS

def format_column(column, spacing):
    # First, get the uniform length
    longest = 0
    for item in column:
        length = len(item)+spacing
        if length > longest:
            longest = length
    # Now format
    return [item+(' '*(longest-len(item))) for item in column]

def format_columns(column_names, key_order, dicts, spacing = 3):
    # Build key column
    columns = [['', '']]
    columns[0].extend(key_order)
    columns[0] = format_column(columns[0], spacing)
    
    # Notes: Python automatically sorts dictionary contents by key, so this will work.
    # TODO: Nevertheless, make this code less hacky
    
    for i, dic in enumerate(dicts):
        new_column = [column_names[i], '-'*len(column_names[i])]
        new_column.extend([str(dic[key]) for key in key_order])
        columns.append(format_column(new_column, spacing))
        
    # Finally, stringify the table
    tbl_str = ''
    for i in range(len(columns[0])): # For each row
        for column in columns: # For each column
            tbl_str += column[i]
        tbl_str += '\n'
    return tbl_str

def validate_dicts(old, new):
    oldkeys = old.keys()
    newkeys = new.keys()
    for key in newkeys:
        if key not in oldkeys:
            return False
    return True

def format_percent(value):
    # value should be a float between 0 and 1
    trimmed = round(value, 3)
    return trimmed * 100



