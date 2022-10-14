#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse


taxo_level_corresp = {
    'k': 'kingdom',
    'p': 'phylum',
    'c': 'class',
    'o': 'order',
    'f': 'family',
    'g': 'genus',
    's': 'species',
    't': 'strains'}
    

def write_taxo_abundance(output_files, level, taxo, abundance):
    if level not in taxo_level_corresp:
        raise ValueError(level + ' is not a know taxonomic level')
    f_n = taxo_level_corresp[level]
    output_files[f_n].write(taxo + '\t')
    output_files[f_n].write(abundance + '\n')


def format_metaphlan2_output(args):
    taxo_levels_abund_f = {
        'kingdom': open(args.kingdom_abundance_file, 'w'),
        'phylum': open(args.phylum_abundance_file, 'w'),
        'class': open(args.class_abundance_file, 'w'),
        'order': open(args.order_abundance_file, 'w'),
        'family': open(args.family_abundance_file, 'w'),
        'genus': open(args.genus_abundance_file, 'w'),
        'species': open(args.species_abundance_file, 'w'),
        'strains': open(args.strains_abundance_file, 'w')
    }

    for taxo_level_f in taxo_levels_abund_f:
        s = taxo_level_f + '\t' + 'abundance\n'
        taxo_levels_abund_f[taxo_level_f].write(s)

    with open(args.metaphlan2_output, 'r') as input_f:
        with open(args.all_taxo_level_abundance_file, 'w') as output_f:
            s = "kingdom\tphylum\tclass\torder\tfamily\t"
            s += "genus\tspecies\tstrains\tabundance\n"
            output_f.write(s)

            levels_number = 8
            for line in input_f.readlines():
                if line.startswith("#"):
                    continue

                split_line = line[:-1].split('\t')
                all_taxo = split_line[0]
                abundance = split_line[1]

                split_taxo = all_taxo.split('|')
                for level in split_taxo:
                    taxo = level.split('__')[1]
                    taxo = taxo.replace("_", " ")
                    output_f.write(taxo + '\t')

                for i in range(len(split_taxo), levels_number):
                    output_f.write('\t')

                output_f.write(abundance + "\n")

                last_taxo_level = split_taxo[-1].split('__')
                taxo = last_taxo_level[1].replace("_", " ")
                level = last_taxo_level[0]
                write_taxo_abundance(
                    taxo_levels_abund_f,
                    level,
                    taxo,
                    abundance)

    for taxo_level_f in taxo_levels_abund_f:
        taxo_levels_abund_f[taxo_level_f].close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--metaphlan2_output', required=True)
    parser.add_argument('--all_taxo_level_abundance_file', required=True)
    parser.add_argument('--kingdom_abundance_file', required=True)
    parser.add_argument('--phylum_abundance_file', required=True)
    parser.add_argument('--class_abundance_file', required=True)
    parser.add_argument('--order_abundance_file', required=True)
    parser.add_argument('--family_abundance_file', required=True)
    parser.add_argument('--genus_abundance_file', required=True)
    parser.add_argument('--species_abundance_file', required=True)
    parser.add_argument('--strains_abundance_file', required=True)
    args = parser.parse_args()

    format_metaphlan2_output(args)
